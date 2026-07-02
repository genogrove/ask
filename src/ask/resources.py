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
    version="0.6.2",
    git_rev="56602a4aef8059cd4bf31d34ac80e5a868c0a122",
    git_tag="v0.6.2",
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
    """A pinned genomic dataset in the curated catalog."""

    name: str
    url: str
    sha256: str
    description: str = ""


# Curated dataset catalog. Each entry pins an *immutable* release (an explicit
# version, never a "latest"/"current" symlink) by URL + sha256. Only names listed
# here are ever fetched.
RESOURCES: dict[str, Resource] = {
    "gencode.human": Resource(
        name="gencode.human",
        url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz",
        sha256="2aaf245c91ed00e80920953add6cfaffcccc876dc0aceeb6ca0c86d15875899a",
        description="GENCODE v50 comprehensive gene annotation, GRCh38 (GFF3, gzip, 1-based).",
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
    dest = _CACHE / res.sha256 / Path(urlsplit(res.url).path).name
    if dest.exists():
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    tmp = tempfile.NamedTemporaryFile(dir=dest.parent, delete=False)
    tmp_path = Path(tmp.name)
    try:
        with tmp, urllib.request.urlopen(res.url) as resp:  # noqa: S310 — pinned catalog URL
            for chunk in iter(lambda: resp.read(1 << 20), b""):
                digest.update(chunk)
                tmp.write(chunk)
        if digest.hexdigest() != res.sha256:
            raise ValueError(
                f"checksum mismatch for {name!r}: expected {res.sha256}, "
                f"got {digest.hexdigest()}"
            )
        tmp_path.replace(dest)  # atomic within the same directory
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
    """A bgzip-compressed, coordinate-sorted, tabix-indexed copy of ``name``'s GFF.

    Built once from the plain-gzip download (``resolve``) and cached next to it;
    region reads (``pg.GffReader(path, region=...)``) require this. Needs htslib's
    ``bgzip`` and ``tabix`` on PATH. Raises ``RuntimeError`` if they're missing.
    """
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
    """True if ``name``'s bgzip+tabix index already exists (so it won't rebuild)."""
    return (_CACHE / RESOURCES[name].sha256 / "indexed.gff3.gz.tbi").exists()


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
