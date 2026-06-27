# SPDX-License-Identifier: GPL-3.0-or-later
"""GFF -> GffGrove load path. Runs only where pygenogrove is installed (CI);
skipped in the bare skeleton env."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from ask.gff import load_gff

GFF3 = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG1\n"
    "chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tID=exon1;Parent=ENSG1\n"
    "chr2\tHAVANA\tgene\t5000\t6000\t.\t-\t.\tID=ENSG2\n"
)


def _write(tmp_path):
    p = tmp_path / "mini.gff3"
    p.write_text(GFF3)
    return p


def test_type_filter_keeps_only_genes(tmp_path) -> None:
    p = _write(tmp_path)
    assert load_gff(p, types={"gene"}).size() == 2  # the exon is filtered out
    assert load_gff(p).size() == 3


def test_intersect_finds_overlapping_gene(tmp_path) -> None:
    genes = load_gff(_write(tmp_path), types={"gene"})
    inside = pg.GenomicCoordinate("*", 1400, 1600)  # within the chr1 gene
    hits = list(genes.intersect(inside, "chr1"))
    assert len(hits) == 1
    assert hits[0].data["type"] == "gene"  # JSON payload on the universal Grove
    assert hits[0].data["id"] == "ENSG1"
    outside = pg.GenomicCoordinate("*", 3000, 3100)  # past the chr1 gene
    assert len(list(genes.intersect(outside, "chr1"))) == 0


def test_seqid_filter(tmp_path) -> None:
    genes = load_gff(_write(tmp_path), types={"gene"}, seqids={"chr2"})
    assert genes.size() == 1
