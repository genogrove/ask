# SPDX-License-Identifier: GPL-3.0-or-later
"""Load a GFF/GTF file into a queryable universal pygenogrove ``Grove``.

The canonical GENCODE -> Grove transform. Pair with ``resources.resolve``::

    from ask import gff, resources
    g = gff.load_gff(resources.resolve("gencode.human"), region=("chr7", 55_000_000, 55_300_000))

Targets the universal ``pg.Grove`` (JSON payloads + labelled edges), not the
typed ``GffGrove``, so the same grove can later hold non-GFF data (BED enhancers,
a link table) and regulatory edges alongside the GFF structure. The universal
grove takes an explicit ``GenomicCoordinate``, so we do the GFF 1-based-inclusive
-> 0-based-closed conversion here ([start-1, end-1]) — the exact rule pygenogrove's
entry-deriving insert uses internally.

GFF3's gene -> transcript -> exon hierarchy (column-9 ``ID`` / ``Parent``) is
reconstructed as **directed edges** with three relations:

* ``{"rel": "contains"}`` — fully-enumerable structural children: gene -> each
  transcript, transcript -> each CDS/UTR. ``get_neighbors`` gives them all.
* ``{"rel": "first_exon"}`` — transcript -> its 5' exon, the splice-path entry
  (NOT generic containment: it reaches one exon, not all of them).
* ``{"rel": "next"}`` — a transcript's exons chained 5'->3' (strand-aware) from
  that first exon; junctions and introns (the gaps) derive from this chain.

So enumerate an isoform's exons by ``first_exon`` then walking ``next``; the
labels keep containment, splice-order, and later regulatory edges distinguishable.
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

    Each feature becomes a key with a JSON payload
    ``{"type", "id", "name", "biotype"}`` (``id`` = column-9 ``ID``, ``name`` =
    ``gene_name``, ``biotype`` = ``gene_type``; ``None`` when absent). The
    GFF3 ``ID``/``Parent`` hierarchy becomes directed edges: ``{"rel": "contains"}``
    for fully-enumerable children, ``{"rel": "first_exon"}`` for transcript -> 5'
    exon, and ``{"rel": "next"}`` chaining exons 5'->3' (see the module docstring).

    Loads only the slice you ask for — GENCODE has millions of features, so
    filter (``None`` = no filter on that axis; avoid all-None on full GENCODE):

    * ``types`` — keep only these feature types, e.g. ``{"gene"}``. NOTE: keep a
      parent's type too, or its children's containment edges won't form.
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
    by_id: dict[str, object] = {}  # GFF3 ID -> Key, for resolving Parent references
    pending: list[tuple[object, list[str]]] = []  # (child_key, parent_ids)
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
        fid = e.get_attribute("ID")
        # ponytail: fixed payload; thread a payload_fn(entry) if callers need more attrs.
        key = g.insert(e.seqid, coord, {
            "type": e.type, "id": fid,
            "name": e.get_gene_name(), "biotype": e.get_gene_biotype(),
        })
        if fid is not None and fid not in by_id:
            by_id[fid] = key  # gene/transcript IDs are unique; first wins on shared-ID leaves
        parent = e.get_attribute("Parent")
        if parent is not None:
            pending.append((key, parent.split(",")))  # GFF3 allows multiple Parents
    # Resolve edges once every key exists — GFF3 doesn't guarantee parent-before-child.
    exons_by_parent: dict[str, list] = {}  # transcript id -> its exon keys
    for child_key, parent_ids in pending:
        is_exon = child_key.data["type"] == "exon"
        for pid in parent_ids:
            if is_exon:
                # Deferred: reached via transcript -> first exon + the splice chain
                # below, not a containment edge to every exon.
                exons_by_parent.setdefault(pid, []).append(child_key)
                continue
            parent_key = by_id.get(pid)
            if parent_key is not None:  # skip Parents filtered out or absent (dangling edge)
                g.add_edge(parent_key, child_key, {"rel": "contains"})

    # Splice chain: order each transcript's exons 5'->3' ('+' ascending, '-'
    # descending), link the transcript to the FIRST exon (enough to reach the
    # isoform — walk {"rel": "next"} for the rest), then chain them. Only exons
    # are chained; chaining same-type siblings generally would wrongly sequence a
    # gene's alternative transcripts.
    for pid, exons in exons_by_parent.items():
        exons.sort(key=lambda k: k.value.start, reverse=(exons[0].value.strand == "-"))
        parent_key = by_id.get(pid)
        if parent_key is not None:
            g.add_edge(parent_key, exons[0], {"rel": "first_exon"})
        for a, b in zip(exons, exons[1:]):
            g.add_edge(a, b, {"rel": "next"})
    return g
