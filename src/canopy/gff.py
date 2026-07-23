# SPDX-License-Identifier: GPL-3.0-or-later
"""Load a GFF/GTF file into a queryable universal pygenogrove ``Grove``.

The canonical GENCODE -> Grove transform. Pair with ``resources.resolve``::

    from canopy import gff, resources
    g = gff.load_gff(resources.resolve("gencode.human"), region=("chr7", 55_000_000, 55_300_000))

Targets the universal ``pg.Grove`` (JSON payloads + labelled edges), not the
typed ``GffGrove``, so the same grove can later hold non-GFF data (BED enhancers,
a link table) and regulatory edges alongside the GFF structure. The universal
grove takes an explicit ``GenomicCoordinate``, so we do the GFF 1-based-inclusive
-> 0-based-closed conversion here ([start-1, end-1]) — the exact rule pygenogrove's
entry-deriving insert uses internally.

GFF3's gene -> transcript -> exon hierarchy (column-9 ``ID`` / ``Parent``) is
reconstructed as **directed edges** with three relations:

* ``{"rel": "contains"}`` — fully-enumerable structural children, e.g. gene ->
  each transcript. ``get_neighbors`` gives them all.
* ``{"rel": "first_exon"}`` — transcript -> its 5' exon, the splice-path entry
  (NOT generic containment: it reaches one exon, not all of them).
* ``{"rel": "next"}`` — a transcript's exons chained 5'->3' (strand-aware) from
  that first exon; junctions and introns (the gaps) derive from this chain.

**Coding structure is folded into exons, not stored as separate nodes.** CDS is a
single contiguous span per transcript (start -> stop codon), sliced across exons;
5'/3' UTR are the exonic parts outside it. So each exon payload carries
``cds = [start, end] | None`` (its coding sub-range; ``None`` = fully UTR), the
transcript carries ``cds_start`` / ``cds_end``, and UTRs are *derived* (exon minus
``cds``, 5'/3' by strand) — never stored. CDS / UTR / codon GFF features become
this annotation; they are never inserted as keys.

So enumerate an isoform's exons by ``first_exon`` then walking ``next``; the
labels keep containment, splice-order, and later regulatory edges distinguishable.
"""

from __future__ import annotations

from collections.abc import Container
from pathlib import Path

# Translation boundaries are captured as the per-transcript CDS span + per-exon
# coding sub-range, so these GFF feature types are consumed for annotation (CDS)
# or dropped (UTR/codon — derived), never stored as keys.
_DROP_TYPES = frozenset({
    "five_prime_UTR", "three_prime_UTR", "UTR", "start_codon", "stop_codon",
    "Selenocysteine",
})


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
    ``gene_name``, ``biotype`` = ``gene_type``; ``None`` when absent). Transcripts
    also carry ``cds_start`` / ``cds_end`` and exons carry ``cds`` (see the module
    docstring). The ``ID``/``Parent`` hierarchy becomes ``contains`` / ``first_exon``
    / ``next`` edges.

    Loads only the slice you ask for — GENCODE has millions of features, so
    filter (``None`` = no filter on that axis; avoid all-None on full GENCODE):

    * ``types`` — keep only these feature types, e.g. ``{"gene"}``. NOTE: keep a
      parent's type too, or its children's edges won't form. CDS is always read
      for exon annotation regardless of this filter; UTR/codon are always dropped.
    * ``seqids`` — keep only these chromosomes.
    * ``region`` — ``(seqid, start, end)``, 0-based closed (same convention as
      ``GenomicCoordinate``); keep only features overlapping that window. Finer
      than ``seqids`` — load a single locus.

    The file is streamed once into memory, then assembled (CDS spans must be known
    before exons are inserted, and there's no payload-update API). The filters
    bound what's buffered and the grove's size, not the read.
    """
    return _assemble(*_parse(
        path, types=types, seqids=seqids, region=region,
        skip_invalid_lines=skip_invalid_lines,
    ))


def _parse(path, *, types, seqids, region, skip_invalid_lines):
    """Stream the file once: buffer node-feature tuples + each transcript's CDS span."""
    import pygenogrove as pg

    rseqid, rstart, rend = region if region is not None else (None, None, None)

    def passes(seqid, start, end):
        if seqids is not None and seqid not in seqids:
            return False
        if region is not None and (seqid != rseqid or start > rend or end < rstart):
            return False
        return True

    feats = []  # (seqid, start, end, strand, type, id, name, biotype, parent_ids)
    cds_span: dict[str, tuple[int, int]] = {}  # transcript id -> (min, max) 0-based closed
    for e in pg.GffReader(str(path), skip_invalid_lines=skip_invalid_lines):
        start, end = e.start - 1, e.end - 1  # GFF 1-based inclusive -> 0-based closed
        if not passes(e.seqid, start, end):
            continue
        if e.type == "CDS":  # fold into exons; never a node
            parent = e.get_attribute("Parent") or ""
            for pid in parent.split(","):
                if pid:
                    lo, hi = cds_span.get(pid, (start, end))
                    cds_span[pid] = (min(lo, start), max(hi, end))
            continue
        if e.type in _DROP_TYPES:  # UTR/codon: derived, never stored
            continue
        if types is not None and e.type not in types:
            continue
        parent = e.get_attribute("Parent")
        feats.append((
            e.seqid, start, end, e.strand, e.type, e.get_attribute("ID"),
            e.get_gene_name(), e.get_gene_biotype(),
            parent.split(",") if parent else [],
        ))
    return feats, cds_span


def _assemble(feats, cds_span):
    """Build a universal ``pg.Grove`` from parsed features: insert keys (with coding
    annotation baked in), then the contains / first_exon / next edges."""
    import pygenogrove as pg

    g = pg.Grove(order=100)
    by_id: dict[str, object] = {}  # GFF3 ID -> Key, for resolving Parent references
    pending: list[tuple[str, object]] = []  # (parent_id, child_key) for non-exon children
    exons_by_parent: dict[str, list] = {}  # transcript id -> [(start, strand, Key)]
    for seqid, start, end, strand, ftype, fid, name, biotype, parent_ids in feats:
        payload = {"type": ftype, "id": fid, "name": name, "biotype": biotype}
        if ftype == "transcript":
            lo_hi = cds_span.get(fid)
            payload["cds_start"], payload["cds_end"] = lo_hi if lo_hi else (None, None)
        elif ftype == "exon":
            payload["cds"] = _exon_cds(start, end, cds_span, parent_ids)
        key = g.insert(seqid, pg.GenomicCoordinate(strand, start, end), payload)
        if fid is not None and fid not in by_id:
            by_id[fid] = key  # gene/transcript IDs are unique; first wins on shared-ID leaves
        if ftype == "exon":
            for pid in parent_ids:
                exons_by_parent.setdefault(pid, []).append((start, strand, key))
        else:
            for pid in parent_ids:
                pending.append((pid, key))

    # Resolve edges — GFF3 doesn't guarantee parent-before-child, but every key now exists.
    for pid, child_key in pending:
        parent_key = by_id.get(pid)
        if parent_key is not None:  # skip Parents filtered out or absent (dangling edge)
            g.add_edge(parent_key, child_key, {"rel": "contains"})

    # Splice chain: order each transcript's exons 5'->3' ('+' ascending, '-'
    # descending), link the transcript to the FIRST exon, then chain the rest.
    for pid, exons in exons_by_parent.items():
        exons.sort(key=lambda se: se[0], reverse=(exons[0][1] == "-"))
        parent_key = by_id.get(pid)
        if parent_key is not None:
            g.add_edge(parent_key, exons[0][2], {"rel": "first_exon"})
        for (_, _, a), (_, _, b) in zip(exons, exons[1:]):
            g.add_edge(a, b, {"rel": "next"})
    return g


def write_sharded_groves(path, out_dir, *, types=None, skip_invalid_lines=False):
    """Parse ``path`` once and serialize a whole-genome ``_all.gg`` plus one
    ``<seqid>.gg`` per chromosome into ``out_dir``; return the seqids written.

    Per-chromosome shards let a query deserialize only the chromosome(s) it needs;
    ``_all.gg`` is the whole-genome grove for genome-wide or cross-chromosome
    queries. The GFF hierarchy edges are all within-chromosome, so sharding by
    chromosome splits no edge.
    """
    out_dir = Path(out_dir)
    feats, cds_span = _parse(
        path, types=types, seqids=None, region=None, skip_invalid_lines=skip_invalid_lines
    )
    _assemble(feats, cds_span).serialize(str(out_dir / "_all.gg"))
    seqids = sorted({f[0] for f in feats})
    for sid in seqids:
        _assemble([f for f in feats if f[0] == sid], cds_span).serialize(str(out_dir / f"{sid}.gg"))
    return seqids


def build_grove(gff_path, region=""):
    """Build the modelled Grove from a bgzip+tabix GFF, reading only ``region``.

    ``region`` is a tabix string (1-based inclusive), e.g. ``"chr7:55000000-55300000"``;
    ``""`` reads the whole file. Returns the same universal Grove as ``load_gff``
    (gene/transcript/exon keys, contains/first_exon/next edges, cds folded onto
    exons). Only features **overlapping** ``region`` are loaded — pick a region that
    covers the features your query needs (a point for "what overlaps here", a gene's
    span for its full structure).

    Self-contained (only ``pygenogrove``) on purpose: its source is injected into the
    sandbox so generated code can call it without importing ``ask``.
    """
    import pygenogrove as pg

    drop = {"five_prime_UTR", "three_prime_UTR", "UTR", "start_codon", "stop_codon", "Selenocysteine"}
    feats, cds = [], {}
    for e in pg.GffReader(str(gff_path), region=region):
        s, en = e.start - 1, e.end - 1  # GFF 1-based inclusive -> 0-based closed
        if e.type == "CDS":
            for pid in (e.get_attribute("Parent") or "").split(","):
                if pid:
                    lo, hi = cds.get(pid, (s, en))
                    cds[pid] = (min(lo, s), max(hi, en))
            continue
        if e.type in drop:
            continue
        parent = e.get_attribute("Parent")
        feats.append((e.seqid, s, en, e.strand, e.type, e.get_attribute("ID"),
                      e.get_gene_name(), e.get_gene_biotype(), parent.split(",") if parent else []))

    g = pg.Grove(order=100)
    by_id, pending, exons_by_parent = {}, [], {}
    for seqid, s, en, strand, ftype, fid, name, biotype, pids in feats:
        pl = {"type": ftype, "id": fid, "name": name, "biotype": biotype}
        if ftype == "transcript":
            sp = cds.get(fid)
            pl["cds_start"], pl["cds_end"] = sp if sp else (None, None)
        elif ftype == "exon":
            pl["cds"] = next(([max(s, cds[p][0]), min(en, cds[p][1])]
                              for p in pids if p in cds and s <= cds[p][1] and en >= cds[p][0]), None)
        k = g.insert(seqid, pg.GenomicCoordinate(strand, s, en), pl)
        if fid and fid not in by_id:
            by_id[fid] = k
        if ftype == "exon":
            for p in pids:
                exons_by_parent.setdefault(p, []).append((s, strand, k))
        else:
            for p in pids:
                pending.append((p, k))
    for p, ck in pending:
        pk = by_id.get(p)
        if pk is not None:
            g.add_edge(pk, ck, {"rel": "contains"})
    for p, ex in exons_by_parent.items():
        ex.sort(key=lambda t: t[0], reverse=(ex[0][1] == "-"))
        pk = by_id.get(p)
        if pk is not None:
            g.add_edge(pk, ex[0][2], {"rel": "first_exon"})
        for (_, _, a), (_, _, b) in zip(ex, ex[1:]):
            g.add_edge(a, b, {"rel": "next"})
    return g


def _exon_cds(start: int, end: int, cds_span: dict, parent_ids: list) -> list | None:
    """The exon's coding sub-range = its intersection with the transcript's CDS
    span (0-based closed), or ``None`` if the exon is entirely UTR / non-coding."""
    for pid in parent_ids:
        span = cds_span.get(pid)
        if span and start <= span[1] and end >= span[0]:
            return [max(start, span[0]), min(end, span[1])]
    return None
