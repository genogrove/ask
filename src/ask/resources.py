# SPDX-License-Identifier: GPL-3.0-or-later
"""Curated resource catalog — the Level 2 reproducibility layer.

A run is reproducible when the question, the resolved datasets, and the library
builds are all pinned. This module is the single source of truth for those pins:

* **Datasets** — named genomic resources with a pinned URL and checksum, resolved
  to a local path on demand (with checksum verification).
* **Builds** — the exact ``pygenogrove`` / ``genogrove`` versions a run was made
  against, recorded so results can be regenerated.

Open-web resource discovery (Level 3) is intentionally out of scope. ``resolve``
takes a curated *name*, never a URL — so the only data ever fetched is what
``RESOURCES`` explicitly defines.
"""

from __future__ import annotations

import hashlib
import shlex
import shutil
import subprocess
import tempfile
import urllib.request
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlsplit


# --------------------------------------------------------------------------- #
# Builds — the pinned library versions a run is made against.
#
# pygenogrove is resolved from a pinned git commit (see [tool.uv.sources] in
# pyproject.toml). The pin below mirrors that commit so a run can record, and
# verify against, the exact build it was made with.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class BuildPin:
    """An exact, immutable library build a run is reproducible against."""

    name: str
    version: str  # expected package version (``pygenogrove.__version__``)
    git_rev: str  # immutable commit the build is pinned to
    git_tag: str = ""  # human-readable tag at that commit, if any


# Keep in lockstep with [tool.uv.sources] `rev` in pyproject.toml.
PYGENOGROVE = BuildPin(
    name="pygenogrove",
    version="0.7.2",
    git_rev="2584321499cb23f814456236f1c3564a2efa956c",
    git_tag="v0.7.2",
)


def verify_pygenogrove_build() -> str:
    """Check that the installed ``pygenogrove`` matches the pinned build.

    Imports lazily so this module is usable without ``pygenogrove`` installed
    (e.g. in the skeleton test env). Returns the underlying C++ engine version
    (``pygenogrove.__genogrove_version__``) so a run can record it. Raises
    ``RuntimeError`` on version drift from the pin.
    """
    import pygenogrove

    installed = getattr(pygenogrove, "__version__", None)
    if installed != PYGENOGROVE.version:
        raise RuntimeError(
            f"pygenogrove build drift: pinned {PYGENOGROVE.version} "
            f"(rev {PYGENOGROVE.git_rev}), installed {installed!r}. "
            "Run `uv sync` to match the pinned build."
        )
    return str(getattr(pygenogrove, "__genogrove_version__", ""))


def build_manifest() -> dict[str, str]:
    """Provenance record of the build a run was made against.

    Combines the static pin with the engine version observed at runtime, for
    embedding in a run's output so results can be regenerated.
    """
    return {
        "pygenogrove_version": PYGENOGROVE.version,
        "pygenogrove_git_rev": PYGENOGROVE.git_rev,
        "pygenogrove_git_tag": PYGENOGROVE.git_tag,
        "genogrove_engine_version": verify_pygenogrove_build(),
    }


# --------------------------------------------------------------------------- #
# Datasets — pinned genomic resources (URL + checksum), resolved to a local path.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Resource:
    """A pinned genomic dataset in the curated catalog.

    ``url`` + ``sha256`` pin the annotation file. For region access it must be
    **bgzip-compressed and tabix-indexed**; when the ``.tbi`` is hosted too,
    ``index_url`` + ``index_sha256`` pin it and ``resolve``/``indexed_path`` fetch
    the pair (no local indexing). ``filename`` overrides the local name when it
    can't be derived from the URL (e.g. Zenodo's ``…/files/<name>/content``).
    """

    name: str
    url: str
    sha256: str
    description: str = ""
    filename: str = ""
    index_url: str = ""
    index_sha256: str = ""
    # Prebuilt whole-genome grove (.gg). When set, the genome-wide path downloads this
    # instead of building locally (minutes) — GroveView reads it lazily. Rebuild + re-pin
    # if the .gg format or the source annotation changes.
    grove_url: str = ""
    grove_sha256: str = ""


# Curated dataset catalog. Each entry pins an *immutable* release (an explicit
# version, never a "latest"/"current" symlink) by URL + sha256. Only names listed
# here are ever fetched.
RESOURCES: dict[str, Resource] = {
    "gencode.human": Resource(
        name="gencode.human",
        # Coordinate-sorted, bgzip+tabix build of GENCODE v50 (Zenodo 21123308),
        # derived from GENCODE v50 (upstream sha 2aaf245c…875899a) — see the record.
        url="https://zenodo.org/api/records/21123308/files/gencode.v50.annotation.sorted.gff3.gz/content",
        sha256="2a87d3a39f9e3be6f0c49359724223ba5e0a094f2fc059b2655635888bb223f5",
        filename="gencode.v50.annotation.sorted.gff3.gz",
        index_url="https://zenodo.org/api/records/21123308/files/gencode.v50.annotation.sorted.gff3.gz.tbi/content",
        index_sha256="52020642c93f01c24488d98b446d705a655d31ea39339fad36cced3b9cc9480a",
        # Prebuilt grove: gene/transcript/exon + contains/first_exon/next edges, CDS folded
        # onto exons. pygenogrove v0.7.2, format 0.2. ~90 MB vs the 237 MB tabix GFF.
        grove_url="https://zenodo.org/api/records/21459419/files/gencode.v50.annotation.grove-fmt0.2.gg/content",
        grove_sha256="df0fca51476d974369db97159a7d4431bfa870d275eb159f4cb21466d3d1a47e",
        description="GENCODE v50 comprehensive gene annotation, GRCh38 (GFF3, sorted + bgzip + tabix).",
    ),
}


# Content-addressed cache: <CACHE>/<sha256>/<filename>. A file only lands here
# after its checksum is verified, so a cache hit needs no re-verification.
_CACHE = Path.home() / ".cache" / "genogrove-ask"


def resolve(name: str) -> Path:
    """Resolve a curated resource to a verified local file path.

    ``name`` is a catalog key (``KeyError`` if not curated) — never a URL, so the
    only data ever fetched is what ``RESOURCES`` defines. On a cache miss, streams
    the pinned URL to a temp file while hashing, and only commits it to the cache
    if the sha256 matches. A mismatch is a hard failure: the partial download is
    discarded and nothing is cached.
    """
    res = RESOURCES[name]
    fname = res.filename or Path(urlsplit(res.url).path).name
    return _download(res.url, res.sha256, _CACHE / res.sha256 / fname)


def _download(url: str, sha256: str, dest: Path) -> Path:
    """Stream ``url`` to ``dest`` (cache hit = no-op), verifying its sha256.

    A mismatch is a hard failure: the partial download is discarded, nothing is
    committed. The commit is an atomic rename within ``dest``'s directory. An empty
    ``sha256`` skips verification — for resolve-on-demand files (ENCODE-rE2G) that are
    not pinned until a run freezes them; curated ``RESOURCES`` always pass a checksum.
    """
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    tmp = tempfile.NamedTemporaryFile(dir=dest.parent, delete=False)
    tmp_path = Path(tmp.name)
    try:
        with tmp, urllib.request.urlopen(url) as resp:  # noqa: S310 — pinned catalog URL
            for chunk in iter(lambda: resp.read(1 << 20), b""):
                digest.update(chunk)
                tmp.write(chunk)
        if sha256 and digest.hexdigest() != sha256:
            raise ValueError(
                f"checksum mismatch for {url!r}: expected {sha256}, got {digest.hexdigest()}"
            )
        tmp_path.replace(dest)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise
    return dest


def data_roots(names: Iterable[str]) -> list[str]:
    """Resolve ``names`` to local file paths for the sandbox's read-only roots."""
    return [str(resolve(n)) for n in names]


# --------------------------------------------------------------------------- #
# Region access — bgzip + tabix so a query reads only its locus (see ask.gff.
# build_grove). GENCODE ships plain gzip, so we recompress + index once.
# --------------------------------------------------------------------------- #


def indexed_path(name: str) -> Path:
    """A bgzip-compressed, coordinate-sorted, tabix-indexed GFF for ``name``.

    If the resource pins a hosted index (``index_url``), download the annotation +
    its ``.tbi`` (no local work). Otherwise fall back to building the index locally
    from the plain-gzip download with htslib's ``bgzip``/``tabix`` — a one-time
    ~minutes step; ``RuntimeError`` if those tools are missing.

    Region reads (``pg.GffReader(path, region=...)``) need the ``.tbi`` next to the
    returned ``.gff3.gz``; both paths are placed accordingly.
    """
    res = RESOURCES[name]
    if res.index_url:  # hosted pair — download, don't build
        gff = resolve(name)  # the sorted-bgzip annotation
        _download(res.index_url, res.index_sha256, gff.with_name(gff.name + ".tbi"))
        return gff

    src = resolve(name)  # plain-gzip download
    out = src.with_name("indexed.gff3.gz")
    tbi = out.with_name(out.name + ".tbi")
    if out.exists() and tbi.exists():
        return out
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"{tool!r} not found — install htslib for region access "
                "(e.g. `brew install htslib` / `apt install tabix`)"
            )
    tmp = out.with_name("indexed.tmp.gff3.gz")
    q = shlex.quote(str(src))
    # Header ('#') lines first, then data sorted by (seqid, start) as tabix requires.
    pipeline = (
        f"{{ gzip -dc {q} | grep '^#' ; gzip -dc {q} | grep -v '^#' | sort -k1,1 -k4,4n ; }} "
        f"| bgzip -c > {shlex.quote(str(tmp))}"
    )
    subprocess.run(pipeline, shell=True, check=True)  # noqa: S602 — our own quoted paths
    subprocess.run(["tabix", "-p", "gff", str(tmp)], check=True)
    Path(str(tmp) + ".tbi").replace(tbi)
    tmp.replace(out)  # commit the index (both parts now in place)
    return out


def is_indexed(name: str) -> bool:
    """True if ``name``'s indexed GFF + `.tbi` are already local (no download/build)."""
    res = RESOURCES[name]
    if res.index_url:  # hosted pair
        fname = res.filename or Path(urlsplit(res.url).path).name
        gff = _CACHE / res.sha256 / fname
        return gff.exists() and gff.with_name(gff.name + ".tbi").exists()
    return (_CACHE / res.sha256 / "indexed.gff3.gz.tbi").exists()  # local build


def _all_grove_gg(name: str) -> Path:
    """Path to the lazily-built whole-genome `.gg` (may not exist yet)."""
    return _CACHE / "groves" / f"{RESOURCES[name].sha256}.{_GROVE_SCHEMA}" / "_all.gg"


def ensure_all_grove(name: str) -> Path:
    """Cache the whole-genome grove (`.gg`) if absent, returning its path.

    Prefers the **pinned prebuilt grove** (``grove_url``): a ~90 MB sha-verified download
    (seconds) instead of a local build. Falls back to building from the annotation
    (``build_grove(region="")`` → serialize) — minutes, only if no grove is pinned. Either
    way it's cached; located queries never trigger this.
    """
    gg = _all_grove_gg(name)
    if gg.exists():
        return gg
    res = RESOURCES[name]
    gg.parent.mkdir(parents=True, exist_ok=True)
    if res.grove_url:  # download the pinned .gg (fast, reproducible)
        return _download(res.grove_url, res.grove_sha256, gg)
    from ask.gff import build_grove  # local fallback — no hosted grove pinned

    tmp = gg.with_name(gg.name + ".tmp")
    build_grove(indexed_path(name), region="").serialize(str(tmp))
    tmp.replace(gg)
    return gg


def grove_view(name: str):
    """Open the whole-genome grove as a lazy ``GroveView`` (downloads the pinned `.gg`
    on first use, else builds it). Serves both located and genome-wide queries — pages in
    only the blocks a query touches; no whole-grove load, no per-query rebuild."""
    import pygenogrove as pg

    return pg.GroveView.open(str(ensure_all_grove(name)))


# Bump when ``ask.gff``'s grove model changes, so a stale `.gg` (valid pygenogrove
# but built from an older schema) is rebuilt rather than silently served.
_GROVE_SCHEMA = "1"


def _grove_dir(name: str) -> Path:
    return _CACHE / "groves" / f"{RESOURCES[name].sha256}.{_GROVE_SCHEMA}"


def grove_index(name: str) -> tuple[dict[str, str], str]:
    """Resolve ``name`` to its sharded grove index, building + caching once.

    Returns ``({seqid: shard_path}, all_path)``: one serialized `.gg` per
    chromosome plus a whole-genome ``_all.gg``. A query deserializes only the
    shard(s) for the chromosome(s) it touches (fast, low-memory); ``_all`` is the
    whole-genome grove for genome-wide or cross-chromosome queries. Built in one
    streaming pass on first use (``ask.gff.write_sharded_groves``) and cached under
    ``<cache>/groves/<sha>.<schema>/``; bump ``_GROVE_SCHEMA`` for model changes.
    """
    from ask.gff import write_sharded_groves

    d = _grove_dir(name)
    if not (d / "_all.gg").exists():
        src = resolve(name)  # downloads + sha256-verifies the source once
        tmp = d.with_name(d.name + ".tmp")
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True, exist_ok=True)
        write_sharded_groves(src, tmp, types={"gene", "transcript", "exon"})
        shutil.rmtree(d, ignore_errors=True)
        tmp.replace(d)  # swap the finished index in atomically
    shards = {p.stem: str(p) for p in d.glob("*.gg") if p.name != "_all.gg"}
    return shards, str(d / "_all.gg")


def grove_path(name: str) -> Path:
    """The whole-genome `.gg` for ``name`` (the ``_all`` grove); builds the index if absent."""
    return Path(grove_index(name)[1])


def is_grove_cached(name: str) -> bool:
    """True if ``name``'s grove index is already built (so resolving it won't rebuild)."""
    return (_grove_dir(name) / "_all.gg").exists()


def load_grove(name: str):
    """Deserialize the whole-genome grove for ``name`` (cached). Self-heals once
    if the cached index won't deserialize."""
    import pygenogrove as pg

    try:
        return pg.Grove.deserialize(str(grove_path(name)))
    except Exception:
        shutil.rmtree(_grove_dir(name), ignore_errors=True)  # nuke the index -> rebuild
        return pg.Grove.deserialize(str(grove_path(name)))


# --------------------------------------------------------------------------- #
# ENCODE-rE2G enhancer→gene predictions — two lazy axes (see the design memo):
#   1. biosample: the catalog is metadata; only a requested biosample's edge BED is
#      ever fetched (1 of ~1,460), the rest cost a catalog row.
#   2. region: a fetched BED is bgzip+tabix-indexed once, then a query reads only the
#      rows overlapping its locus — never the whole ~90k-edge file.
# So a query's footprint = (biosamples it names) × (rows in the region it names).
# --------------------------------------------------------------------------- #

_ENCODE = "https://www.encodeproject.org"

# The rE2G thresholded element-gene-links BED has ~56 columns (mostly model-internal
# `.Feature` inputs); we keep only these, selected BY HEADER NAME — the file's own
# `#`-prefixed header line — not by position, because `Score` is the LAST column and the
# tail count/order varies by rE2G version. `class` is a positional label (col 5); `Score`
# (col 56) is the model's calibrated enhancer→gene confidence. Everything else is either
# derivable from the grove (distance, class) or a model input we don't store.
_RE2G_KEEP = (
    "chr", "start", "end", "class", "TargetGene", "TargetGeneEnsemblID",
    "TargetGeneTSS", "isSelfPromoter", "Score",
)

# Catalog of all ENCODE-rE2G prediction annotations, generated by
# tools/fetch_re2g_catalog.py (one row per biosample; ships with the package).
_RE2G_CATALOG = Path(__file__).parent / "data" / "encode_re2g_catalog.tsv"


def re2g_catalog() -> list[dict[str, str]]:
    """The ENCODE-rE2G biosample catalog (metadata only — no data fetched).

    Each entry maps a biosample to its annotation ``accession``; the agent scopes which
    biosample(s) a query needs (axis 1). Raises ``FileNotFoundError`` with a pointer to
    the fetcher if the catalog hasn't been generated yet.
    """
    if not _RE2G_CATALOG.exists():
        raise FileNotFoundError(
            f"{_RE2G_CATALOG} missing — run `python tools/fetch_re2g_catalog.py` once "
            "to generate the ENCODE-rE2G catalog."
        )
    import csv

    with _RE2G_CATALOG.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def re2g_accessions(biosample_term: str) -> list[str]:
    """Annotation accessions for a biosample term (case-insensitive substring), e.g.
    ``re2g_accessions("prostate")`` — the agent picks from these before fetching any data."""
    t = biosample_term.lower()
    return [e["accession"] for e in re2g_catalog() if t in e["biosample_term"].lower()]


# Ontology-id prefix -> the query axis it represents (UBERON = anatomy/tissue, CL = cell
# type, EFO = cell line / disease, NTR = ENCODE novel term). Lets the agent filter cohorts
# by axis without any external ontology file.
_ONTOLOGY_AXIS = {"UBERON": "tissue", "CL": "cell type", "CLO": "cell line",
                  "EFO": "cell line", "NTR": "novel term"}


def re2g_cohorts() -> list[dict]:
    """The catalog grouped into **cohorts** — one per biosample (its ontology id), folding
    the replicate accessions together. This is the request→biosample association layer: an
    agent picks a cohort from this list by ``name`` / ``axis`` / ``type`` (grounded — only
    declared cohorts exist), and the cohort's ``accessions`` then drive a lazy per-cohort
    grove build (merge replicates → edge support ``n``). Derived from the catalog, so it
    ships nothing new.

    Each cohort: ``ontology_id``, ``name`` (biosample term), ``axis`` (from the ontology
    prefix), ``type`` (biosample_type), ``n_replicates``, ``accessions``. Sorted by replicate
    count (most-replicated first), then name.
    """
    groups: dict[str, dict] = {}
    for e in re2g_catalog():
        g = groups.setdefault(e["biosample_id"], {
            "ontology_id": e["biosample_id"], "name": e["biosample_term"],
            "axis": _ONTOLOGY_AXIS.get(e["biosample_id"].split(":")[0], "other"),
            "type": e["biosample_type"], "accessions": [],
        })
        g["accessions"].append(e["accession"])
    for g in groups.values():
        g["n_replicates"] = len(g["accessions"])
    return sorted(groups.values(), key=lambda g: (-g["n_replicates"], g["name"]))


def _re2g_edge_href(accession: str) -> str:
    """Resolve an rE2G annotation to its thresholded element-gene-links BED download URL.

    Picks the ENCODE-rE2G (not ABC) default file: ``preferred_default`` +
    ``output_type == 'thresholded element gene links'``. Raises if none is found.
    """
    import json

    url = f"{_ENCODE}/annotations/{accession}/?format=json"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as resp:  # noqa: S310 — fixed ENCODE host
        data = json.load(resp)
    for f in data.get("files", []):
        if (f.get("output_type") == "thresholded element gene links"
                and f.get("preferred_default") and f.get("file_format") == "bed"):
            return _ENCODE + f["href"]
    raise RuntimeError(f"no default thresholded rE2G BED found for {accession}")


def re2g_indexed(accession: str) -> Path:
    """Download + bgzip + tabix-index a biosample's rE2G edge BED once; cached.

    Fetches the whole (~12 MB) BED on first use — the region laziness is on the read
    side (``re2g_edges``), not the download. ``RuntimeError`` if htslib is absent.
    """
    out = _CACHE / "re2g" / f"{accession}.bed.gz"
    tbi = out.with_name(out.name + ".tbi")
    if out.exists() and tbi.exists():
        return out
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"{tool!r} not found — install htslib for rE2G region access "
                "(e.g. `brew install htslib` / `apt install tabix`)"
            )
    raw = _download(_re2g_edge_href(accession), "", out.with_name("raw.bed.gz"))  # unpinned
    tmp = out.with_name("indexed.tmp.bed.gz")
    q = shlex.quote(str(raw))
    # Comment/header ('#') lines first, then data sorted by (chrom, start) for tabix.
    pipeline = (
        f"{{ gzip -dc {q} | grep '^#' ; gzip -dc {q} | grep -v '^#' | sort -k1,1 -k2,2n ; }} "
        f"| bgzip -c > {shlex.quote(str(tmp))}"
    )
    subprocess.run(pipeline, shell=True, check=True)  # noqa: S602 — our own quoted paths
    subprocess.run(["tabix", "-p", "bed", str(tmp)], check=True)
    Path(str(tmp) + ".tbi").replace(tbi)
    tmp.replace(out)  # commit both parts together
    raw.unlink(missing_ok=True)
    return out


def _re2g_header(path) -> list[str]:
    """The rE2G file's column names, from its leading ``#``-prefixed header line."""
    import gzip

    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                return ln[1:].rstrip("\n").split("\t")
    raise RuntimeError(f"{path}: no '#'-prefixed header line to name columns")


def re2g_edges(accession: str, region: str = "") -> list[dict[str, str]]:
    """Enhancer→gene edges for a biosample as dicts holding only the ``_RE2G_KEEP``
    columns, selected **by header name** (robust to the ~56-column tail; ``Score`` is last).

    ``region`` is a tabix string (``"chr7:55000000-55300000"``, 1-based inclusive); only
    rows overlapping it are read (axis 2). Empty ``region`` streams the whole file.
    """
    path = re2g_indexed(accession)
    idx = {name: i for i, name in enumerate(_re2g_header(path))}
    if region:
        out = subprocess.run(["tabix", str(path), region],  # noqa: S603
                             capture_output=True, text=True, check=True).stdout
        lines = out.splitlines()
    else:
        import gzip
        with gzip.open(path, "rt") as fh:
            lines = [ln for ln in fh.read().splitlines() if not ln.startswith("#")]
    return [{k: f[idx[k]] for k in _RE2G_KEEP}
            for f in (ln.split("\t") for ln in lines if ln)]


def _tss_pos(raw: str) -> int:
    """Parse a TSS coordinate from an rE2G ``target_tss`` cell — tolerate a bare int,
    ``chr:pos``, or ``start-end`` (take the first)."""
    return int(raw.split("-")[0].split(":")[-1])


def augment_grove(base_gg, edges):
    """Augment the GENCODE grove **in place** with rE2G enhancer→gene edges.

    Deserialize ``base_gg`` (the built GENCODE `.gg`) into a mutable ``pg.Grove`` —
    every gene/transcript/exon key already present — then for each rE2G row add an
    **enhancer** node and a ``{"rel": "regulates"}`` edge from it onto the *existing*
    GENCODE **gene** key. That cross-index edge over a shared gene key is the point:
    a query goes ``variant ∩ enhancer → get_neighbors_if("regulates") → gene →
    first_exon/next`` in one traversal, no second grove.

    The gene is located by intersecting the grove at the rE2G ``target_tss`` and matching
    the GENCODE gene's ENSG id (versioned, so compared on the base) to ``target_ensembl``.
    rE2G targets absent from GENCODE are counted, not invented. BED is half-open; the
    key is 0-based closed, so ``end`` shifts by one. Returns ``(grove, stats)``.
    """
    import pygenogrove as pg

    g = pg.Grove.deserialize(str(base_gg))
    enh: dict[tuple, object] = {}       # (chrom, start, end) -> enhancer Key (deduped)
    gene_cache: dict[str, object] = {}  # base ENSG id -> gene Key, or None (confirmed miss)

    def find_gene(chrom, tss, base):
        if base in gene_cache:
            return gene_cache[base]
        hit = None
        # "*" is the strand wildcard — genes carry real strands (+/-), so a "." query
        # would match nothing. See pygenogrove test_object_grove (strand is significant).
        for k in g.intersect(pg.GenomicCoordinate("*", tss, tss), chrom):
            d = k.data
            if d.get("type") == "gene" and (d.get("id") or "").split(".")[0] == base:
                hit = k
                break
        gene_cache[base] = hit
        return hit

    linked = missed = self_prom = 0
    for e in edges:
        gene = find_gene(e["chr"], _tss_pos(e["TargetGeneTSS"]), e["TargetGeneEnsemblID"].split(".")[0])
        if gene is None:  # rE2G target not in this GENCODE build — report, don't fabricate
            missed += 1
            continue
        # Self-promoters (element IS the gene's own promoter) are kept — they're
        # self-identifying (class="promoter" + ~0 derivable distance), so let the agent
        # filter them, don't drop links at build time. Tallied for visibility.
        if e["isSelfPromoter"].upper() == "TRUE":
            self_prom += 1
        chrom, es, ee = e["chr"], int(e["start"]), int(e["end"]) - 1  # half-open -> closed
        ek = (chrom, es, ee)
        if ek not in enh:  # one node per element; class (col 5) is the only stored annotation
            enh[ek] = g.insert(chrom, pg.GenomicCoordinate(".", es, ee),  # unstranded element
                               {"type": "enhancer", "class": e["class"]})
        # score (col 56) is the whole payload — distance/class are derivable from the nodes.
        g.add_edge(enh[ek], gene, {"rel": "regulates", "score": float(e["Score"])})
        linked += 1
    return g, {"enhancers": len(enh), "regulates": linked,
               "missed_targets": missed, "self_promoters": self_prom}


def ensure_augmented_grove(base_name: str, accession: str) -> Path:
    """Build + cache the combined grove = ``base_name``'s GENCODE grove augmented with
    biosample ``accession``'s rE2G edges. Returns its `.gg` path (built once, then cached
    beside the base grove). Fetching the edges needs htslib + network; the query path only
    opens the resulting local `.gg` via ``GroveView``.
    """
    gg = _all_grove_gg(base_name).with_name(f"+re2g-{accession}.gg")
    if gg.exists():
        return gg
    base_gg = ensure_all_grove(base_name)  # the pinned GENCODE .gg
    gg.parent.mkdir(parents=True, exist_ok=True)
    tmp = gg.with_name(gg.name + ".tmp")
    grove, _stats = augment_grove(base_gg, re2g_edges(accession, ""))
    grove.serialize(str(tmp))
    tmp.replace(gg)
    return gg
