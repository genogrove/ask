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
