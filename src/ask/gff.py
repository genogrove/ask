# SPDX-License-Identifier: GPL-3.0-or-later
"""Load a GFF/GTF file into a queryable universal pygenogrove ``Grove``.

The canonical GENCODE -> Grove transform. Pair with ``resources.resolve``::

    from ask import gff, resources
    genes = gff.load_gff(resources.resolve("gencode.human"), types={"gene"})

Targets the universal ``pg.Grove`` (JSON payloads), not the typed ``GffGrove``,
so the same grove can later hold non-GFF data (BED enhancers, a link table) and
labelled edges. The universal grove takes an explicit ``GenomicCoordinate``, so
we do the GFF 1-based-inclusive -> 0-based-closed conversion here ([start-1,
end-1]) — the exact rule pygenogrove's entry-deriving insert uses internally.
"""

from __future__ import annotations

from collections.abc import Container
from pathlib import Path


def load_gff(
    path: str | Path,
    *,
    types: Container[str] | None = None,
    seqids: Container[str] | None = None,
    region: tuple[str, int, int] | None = None,
    skip_invalid_lines: bool = False,
):
    """Read ``path`` (plain/gzip/BGZF GFF or GTF) into a universal ``pg.Grove``.

    Each feature becomes a key with a JSON payload ``{"type", "id", "name"}``
    (``id`` = column-9 ``ID``, ``name`` = ``gene_name`` attribute; ``None`` when
    absent). Loads only the slice you ask for — GENCODE has millions of features,
    so filter (``None`` = no filter on that axis; avoid all-None on full GENCODE):

    * ``types`` — keep only these feature types, e.g. ``{"gene"}``.
    * ``seqids`` — keep only these chromosomes.
    * ``region`` — ``(seqid, start, end)``, 0-based closed (same convention as
      ``GenomicCoordinate``); keep only features overlapping that window. Finer
      than ``seqids`` — load a single locus.

    The file is always streamed in full (the gzip isn't tabix-indexed, so there's
    no random access); the filters bound the grove's size, not the read.
    """
    import pygenogrove as pg

    rseqid, rstart, rend = region if region is not None else (None, None, None)
    g = pg.Grove(order=100)
    for e in pg.GffReader(str(path), skip_invalid_lines=skip_invalid_lines):
        if types is not None and e.type not in types:
            continue
        if seqids is not None and e.seqid not in seqids:
            continue
        # GFF 1-based inclusive [start, end] -> 0-based closed [start-1, end-1].
        start, end = e.start - 1, e.end - 1
        if region is not None and (e.seqid != rseqid or start > rend or end < rstart):
            continue
        coord = pg.GenomicCoordinate(e.strand, start, end)
        # ponytail: fixed payload; thread a payload_fn(entry) if callers need more attrs.
        g.insert(e.seqid, coord, {"type": e.type, "id": e.get_attribute("ID"), "name": e.get_gene_name()})
    return g
