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
