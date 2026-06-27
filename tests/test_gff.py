# SPDX-License-Identifier: GPL-3.0-or-later
"""GFF -> universal Grove load path. Runs only where pygenogrove is installed
(CI); skipped in the bare skeleton env."""

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


def test_region_filter(tmp_path) -> None:
    p = _write(tmp_path)  # chr1 gene 1-based [1000,2000] -> 0-based closed [999,1999]
    assert load_gff(p, region=("chr1", 1400, 1600)).size() == 1  # overlaps the gene
    assert load_gff(p, region=("chr1", 2500, 3000)).size() == 0  # window past every chr1 feature
    assert load_gff(p, region=("chr2", 999, 1999)).size() == 0   # right window, wrong chromosome


HIER = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t3000\t.\t+\t.\tID=g1;gene_name=AAA;gene_type=protein_coding\n"
    "chr1\tHAVANA\ttranscript\t1000\t3000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tHAVANA\texon\t2900\t3000\t.\t+\t.\tID=e2;Parent=t1\n"
)


def test_hierarchy_edges(tmp_path) -> None:
    p = tmp_path / "hier.gff3"
    p.write_text(HIER)
    g = load_gff(p)  # all feature types, so the hierarchy is intact

    gene = next(
        k for k in g.intersect(pg.GenomicCoordinate("*", 1500, 1500), "chr1")
        if k.data["type"] == "gene"
    )
    assert gene.data["name"] == "AAA"
    assert gene.data["biotype"] == "protein_coding"

    txs = list(g.get_neighbors(gene))  # gene -> transcript
    assert [t.data["type"] for t in txs] == ["transcript"]
    assert g.get_edges(gene)[0] == {"rel": "contains"}  # labelled containment edge

    # transcript -> first (5') exon only; the rest hang off the splice chain
    entry = list(g.get_neighbors(txs[0]))
    assert [x.data["id"] for x in entry] == ["e1"]
    assert g.get_edges(txs[0])[0] == {"rel": "first_exon"}

    # splice chain: '+' strand, ascending order e1 -> e2
    assert [n.data["id"] for n in g.get_neighbors(entry[0])] == ["e2"]
    assert g.get_edges(entry[0])[0] == {"rel": "next"}

    # no CDS in this fixture -> exons are fully non-coding, transcript span is None
    assert entry[0].data["cds"] is None
    assert txs[0].data["cds_start"] is None


CODING = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t3000\t.\t+\t.\tID=g1;gene_name=AAA;gene_type=protein_coding\n"
    "chr1\tHAVANA\ttranscript\t1000\t3000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tHAVANA\texon\t2800\t3000\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tHAVANA\tCDS\t1100\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tHAVANA\tCDS\t2800\t2900\t.\t+\t2\tID=c2;Parent=t1\n"
    "chr1\tHAVANA\tfive_prime_UTR\t1000\t1099\t.\t+\t.\tID=u1;Parent=t1\n"
)


def test_cds_folded_into_exons(tmp_path) -> None:
    p = tmp_path / "coding.gff3"
    p.write_text(CODING)
    g = load_gff(p)

    everything = list(g.intersect(pg.GenomicCoordinate("*", 0, 5000), "chr1"))
    # CDS and UTR are folded/derived, never inserted as keys
    assert {k.data["type"] for k in everything} == {"gene", "transcript", "exon"}

    tx = next(k for k in everything if k.data["type"] == "transcript")
    assert (tx.data["cds_start"], tx.data["cds_end"]) == (1099, 2899)  # CDS span, 0-based closed

    e1 = next(k for k in everything if k.data["id"] == "e1")  # exon 999..1199, coding from 1099
    assert e1.data["cds"] == [1099, 1199]                     # 999..1098 is 5' UTR (derived)
    e2 = next(k for k in everything if k.data["id"] == "e2")  # exon 2799..2999, coding to 2899
    assert e2.data["cds"] == [2799, 2899]                     # 2900..2999 is 3' UTR (derived)


MINUS = (
    "##gff-version 3\n"
    "chr3\tHAVANA\tgene\t100\t400\t.\t-\t.\tID=mg\n"
    "chr3\tHAVANA\ttranscript\t100\t400\t.\t-\t.\tID=mt;Parent=mg\n"
    "chr3\tHAVANA\texon\t100\t200\t.\t-\t.\tID=lo;Parent=mt\n"
    "chr3\tHAVANA\texon\t300\t400\t.\t-\t.\tID=hi;Parent=mt\n"
)


def test_splice_chain_is_strand_aware(tmp_path) -> None:
    p = tmp_path / "minus.gff3"
    p.write_text(MINUS)
    g = load_gff(p)
    # 5'->3' on the '-' strand runs high coordinate -> low: hi -> lo.
    hi = next(
        k for k in g.intersect(pg.GenomicCoordinate("*", 349, 349), "chr3")
        if k.data["id"] == "hi"
    )
    assert [n.data["id"] for n in g.get_neighbors(hi)] == ["lo"]
