# SPDX-License-Identifier: GPL-3.0-or-later
"""Curated resource registry — the Level 2 reproducibility layer.

A run is reproducible when the question, the resolved datasets, and the library
builds are all pinned. This module is the single source of truth for those pins:

* **Datasets** — named genomic resources with a pinned URL and checksum, resolved
  to a local path on demand (with checksum verification).
* **Builds** — the exact ``pygenogrove`` / ``genogrove`` versions a run was made
  against, recorded so results can be regenerated.

Open-web resource discovery (Level 3) is intentionally out of scope; only entries
present in the curated registry are available to a query.
"""

from __future__ import annotations

from dataclasses import dataclass


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
    version="0.1.0",
    git_rev="70d77ea039cbb4bf16544a4edfd86a689a8720fb",
    git_tag="v0.1.0",
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
    """A pinned genomic dataset in the curated registry."""

    name: str
    url: str
    sha256: str
    description: str = ""


# Curated dataset registry. Populated as hero-query resources are added; each
# entry must carry a pinned URL and checksum.
REGISTRY: dict[str, Resource] = {}


def resolve(name: str) -> Resource:
    """Look up a registry entry by name. Raises ``KeyError`` if not curated."""
    return REGISTRY[name]
