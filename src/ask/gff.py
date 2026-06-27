# SPDX-License-Identifier: GPL-3.0-or-later
"""Load a GFF/GTF file into a queryable pygenogrove ``GffGrove``.

The canonical GENCODE -> Grove transform. Pair with ``resources.resolve``::

    from ask import gff, resources
    genes = gff.load_gff(resources.resolve("gencode.human"), types={"gene"})

Uses ``GffGrove``'s entry-deriving insert so the GFF 1-based coordinates and
strand are converted correctly — don't hand-roll that conversion, it's the
off-by-one trap. A labelled-edge graph needs the universal ``pg.Grove`` (explicit
``GenomicCoordinate`` + JSON payload); that's a separate assembly step.
"""

from __future__ import annotations

from collections.abc import Container
from pathlib import Path


def load_gff(
    path: str | Path,
    *,
    types: Container[str] | None = None,
    seqids: Container[str] | None = None,
    skip_invalid_lines: bool = False,
):
    """Read ``path`` (plain/gzip/BGZF GFF or GTF) into a ``pg.GffGrove``.

    Loads only the slice you ask for — GENCODE has millions of features, so
    filter: ``types`` keeps only those feature types (e.g. ``{"gene"}``),
    ``seqids`` only those chromosomes. ``None`` means no filter on that axis
    (avoid on full GENCODE).
    """
    import pygenogrove as pg

    g = pg.GffGrove(order=100)
    for e in pg.GffReader(str(path), skip_invalid_lines=skip_invalid_lines):
        if types is not None and e.type not in types:
            continue
        if seqids is not None and e.seqid not in seqids:
            continue
        g.insert(e.seqid, e)  # entry-deriving: 1-based->closed + strand from the record
    # ponytail: per-entry insert; switch to insert_bulk(presorted) if GENCODE load time bites.
    return g
