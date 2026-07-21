# SPDX-License-Identifier: GPL-3.0-or-later
"""Augment the GENCODE grove with one ENCODE-rE2G biosample's enhancer→gene edges.

    python tools/build_re2g_grove.py ENCSR864JGD     # an LNCaP clone FGC accession

Resolves the (pinned) GENCODE `.gg`, downloads the biosample's thresholded
element-gene-links BED (once, cached), and writes a *combined* grove: GENCODE
gene/transcript/exon nodes + LNCaP enhancer nodes + `regulates` edges hung on the
existing gene keys. Prints the counts (enhancers / edges / rE2G targets not found in
GENCODE) and the final `.gg` size. This is the local build; the produced `.gg` is what
gets pinned to Zenodo and wired into the query path.

Run `--selftest` to exercise the augmenter on a synthetic 1-gene grove — no network,
no htslib, no GENCODE download.
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ask import resources  # noqa: E402


def _row(**kw):
    """An rE2G edge dict (only the kept columns), overridable by keyword."""
    base = {"chr": "chr8", "start": "0", "end": "0", "class": "intergenic",
            "TargetGene": "GENE", "TargetGeneEnsemblID": "ENSG0", "TargetGeneTSS": "0",
            "isSelfPromoter": "FALSE", "Score": "0.5"}
    base.update(kw)
    return base


def _selftest() -> None:
    """Augment a synthetic grove holding one gene (MYC) and confirm the enhancer +
    regulates edge land on the existing gene key, an off-target is counted, a
    self-promoter is kept + tallied, and it survives a serialize/GroveView round-trip."""
    import pygenogrove as pg

    # A minimal GENCODE-like base grove: one gene node, versioned ENSG id.
    base = pg.Grove(order=100)
    base.insert("chr8", pg.GenomicCoordinate("+", 127735000, 127742000),
                {"type": "gene", "id": "ENSG00000136997.20", "name": "MYC"})

    base_gg = Path(tempfile.mkdtemp()) / "base.gg"
    base.serialize(str(base_gg))

    edges = [
        # two enhancers targeting MYC (unversioned ENSG, TSS inside the gene body)
        _row(start="100", end="200", TargetGene="MYC",
             TargetGeneEnsemblID="ENSG00000136997", TargetGeneTSS="127736000", Score="0.9"),
        _row(start="300", end="400", TargetGene="MYC",
             TargetGeneEnsemblID="ENSG00000136997", TargetGeneTSS="127736000", Score="0.7"),
        # a self-promoter at the TSS -> kept as a regulates edge, and tallied
        _row(start="127735800", end="127736100", TargetGene="MYC",
             TargetGeneEnsemblID="ENSG00000136997", TargetGeneTSS="127736000",
             isSelfPromoter="TRUE"),
        # a target not in the base grove -> counted as missed, not invented
        _row(start="500", end="600", TargetGene="GHOST",
             TargetGeneEnsemblID="ENSG99999999", TargetGeneTSS="127736000"),
    ]
    g, stats = resources.augment_grove(base_gg, edges)
    assert stats == {"enhancers": 3, "regulates": 3,
                     "missed_targets": 1, "self_promoters": 1}, stats

    out = Path(tempfile.mkdtemp()) / "combined.gg"
    g.serialize(str(out))
    gv = pg.GroveView.open(str(out))

    # variant ∩ enhancer -> regulates -> the SAME MYC gene node (shared key).
    # "*" wildcard strand — enhancers are unstranded, genes are +/-.
    hits = list(gv.intersect(pg.GenomicCoordinate("*", 150, 150), "chr8"))
    enh = [h for h in hits if h.data.get("type") == "enhancer"]
    assert len(enh) == 1 and enh[0].data["class"] == "intergenic", [h.data for h in enh]
    assert enh[0].value.end == 199, enh[0].value.end  # half-open [100,200) -> closed [100,199]
    targets = gv.get_neighbors_if(enh[0], lambda m: m.get("rel") == "regulates")
    assert [t.data["name"] for t in targets] == ["MYC"], [t.data for t in targets]
    assert targets[0].data["id"] == "ENSG00000136997.20", targets[0].data  # the GENCODE key
    edges_meta = gv.get_edges(enh[0])
    assert edges_meta[0]["score"] == 0.9, edges_meta  # real Score (col 56), not a feature col
    print("augment_grove selftest OK:", stats)


def main(argv):
    if not argv or argv[0] == "--selftest":
        _selftest()
        return 0
    accession = argv[0]
    print(f"Augmenting GENCODE with rE2G {accession} (GENCODE .gg + edge download + "
          "in-memory deserialize; first run is slow)…", file=sys.stderr)
    base_gg = resources.ensure_all_grove("gencode.human")
    grove, stats = resources.augment_grove(base_gg, resources.re2g_edges(accession, ""))
    out = resources._all_grove_gg("gencode.human").with_name(f"+re2g-{accession}.gg")
    out.parent.mkdir(parents=True, exist_ok=True)
    grove.serialize(str(out))
    print(f"grove: {out}")
    print(f"size:  {out.stat().st_size / 1e6:.1f} MB")
    print(f"stats: {stats['enhancers']} enhancers, {stats['regulates']} regulates edges, "
          f"{stats['missed_targets']} rE2G targets not in GENCODE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
