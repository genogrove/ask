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
    version="0.7.1",
    git_rev="da38261823ea94beda534e313f5c72f061b97618",
    git_tag="v0.7.1",
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
    """Build + cache the whole-genome grove (for genome-wide queries) if absent.

    Slow on first call (reads the whole annotation), instant after. Only invoked
    when a query actually needs the whole genome — located queries never build it.
    """
    from ask.gff import build_grove

    gg = _all_grove_gg(name)
    if not gg.exists():
        gg.parent.mkdir(parents=True, exist_ok=True)
        tmp = gg.with_name(gg.name + ".tmp")
        build_grove(indexed_path(name), region="").serialize(str(tmp))
        tmp.replace(gg)
    return gg


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

# The thresholded element-gene-links BED columns, in file order (see the schema on a
# released ENCODE-rE2G annotation). Used to parse a tabix region slice into dicts.
RE2G_FIELDS = (
    "chrom", "start", "end", "element", "class", "target_gene", "target_ensembl",
    "target_tss", "is_self_promoter", "cell_type", "distance_to_tss", "dnase_prom",
    "contact_3d", "abc_score", "n_candidate_enh_gene", "n_tss_enh_gene",
    "sum_nearby_enh", "ubiquitous_gene", "score",
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


def re2g_edges(accession: str, region: str = "") -> list[dict[str, str]]:
    """Enhancer→gene edges for a biosample, as dicts keyed by ``RE2G_FIELDS``.

    ``region`` is a tabix string (``"chr7:55000000-55300000"``, 1-based inclusive); only
    rows overlapping it are read (axis 2). Empty ``region`` streams the whole file.
    """
    path = re2g_indexed(accession)
    if region:
        out = subprocess.run(["tabix", str(path), region],  # noqa: S603
                             capture_output=True, text=True, check=True).stdout
        lines = out.splitlines()
    else:
        import gzip
        with gzip.open(path, "rt") as fh:
            lines = [ln for ln in fh.read().splitlines() if not ln.startswith("#")]
    return [dict(zip(RE2G_FIELDS, ln.split("\t"))) for ln in lines if ln]
