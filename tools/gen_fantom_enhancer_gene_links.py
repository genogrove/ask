# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate FANTOM5 enhancer -> gene `coactive` edges ourselves (hg38).

The published `enhancer_tss_associations.bed` (Andersson 2014) is hg19 and only
lived on a now-dead host, so we regenerate it from the raw FANTOM5 CAGE matrices,
faithfully to the paper's method but hg38-native and reproducible:

  1. enhancer CAGE signal  (F5.hg38.enhancers.expression.*.matrix, rows = enhancers)
  2. promoter CAGE signal  (hg38 CAGE-peak expression matrix, rows = peaks)
  3. average libraries into facets (optional sample->facet map; strongly recommended
     -- raw ~1800 libraries over-power the correlation and give anticonservative FDR)
  4. for each enhancer, Pearson-correlate against every CAGE peak whose TSS lies
     within +/-500 kb, across facets
  5. Benjamini-Hochberg over ALL tested pairs; keep FDR < cutoff and r > 0
  6. map each surviving promoter peak to a GENCODE gene by TSS proximity, collapse to
     one edge per (enhancer, gene) keeping the strongest correlation

Output is a `--links` edge table: enhancer_id, gene_id, gene_name, rel=coactive, R, FDR
-- the `coactive` edges for the grove. Row-ids in both matrices already encode hg38
coordinates (enhancer `chr:start-end`, peak `chr:start..end,strand`), so no BED needed.

Deliberately gene-target, promoter-mediated -- NOT a substitute for Hi-C `loops_to` or
GTEx `eqtl`: the point of `coactive` is to be an INDEPENDENT modality (co-expression)
so ">=2 concordant modalities -> trust the target" means something. Keep it FANTOM.

Run (heavy: the peak matrix is ~1.5 GB in memory -- run on a real box, not in-session):

    python gen_fantom_enhancer_gene_links.py \
        --enhancer-matrix F5.hg38.enhancers.expression.tpm.matrix.gz \
        --peak-matrix     hg38_CAGE_peaks_expression_tpm.matrix.gz \
        --gencode         gencode.v50.annotation.gff3.gz \
        --facets          library_to_facet.tsv \
        --out             fantom_enhancer_gene.coactive.tsv

Needs numpy + scipy + pandas (a data-prep script, not part of the `ask` runtime).
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

WINDOW = 500_000  # enhancer<->promoter search radius (bp), per Andersson 2014
TSS_WINDOW = 500  # a CAGE peak is a promoter of a gene if within this of its TSS
R_FLOOR = 0.1     # store only pairs with r above this before FDR (memory guard; see note)
FDR = 0.05
ID_RE = r"CNhs\d+"  # FANTOM library id embedded in matrix column headers


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_matrix(path: str, id_regex: str) -> tuple[list[str], np.ndarray, list[str]]:
    """Load a FANTOM expression matrix -> (row_ids, data[float32], library_ids).

    Columns are keyed by the FANTOM library id (``CNhs\\d+``) parsed from each header,
    so the enhancer and peak matrices can be aligned on shared libraries even when
    their full header strings differ. Non-library columns (annotation) are dropped.
    """
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0, low_memory=False)
    rx = re.compile(id_regex)
    keep = {c: m.group(0) for c in df.columns if (m := rx.search(str(c)))}
    df = df[list(keep)].astype("float32")
    df.columns = [keep[c] for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]  # first library column wins on collision
    return list(df.index.astype(str)), df.to_numpy(), list(df.columns)


def facet_average(data: np.ndarray, libs: list[str], mapping: dict[str, str]):
    """Average library columns into facets. Returns (data[n_rows, n_facets], facets)."""
    facets, cols = [], []
    by_facet: dict[str, list[int]] = {}
    for j, lib in enumerate(libs):
        f = mapping.get(lib)
        if f is not None:
            by_facet.setdefault(f, []).append(j)
    for f, idx in by_facet.items():
        facets.append(f)
        cols.append(data[:, idx].mean(axis=1))
    return np.column_stack(cols).astype("float32"), facets


def parse_point(row_id: str) -> tuple[str, int] | None:
    """Representative genomic point from a matrix row-id.

    Enhancer id ``chr1:839741-840250`` -> midpoint. CAGE-peak id
    ``chr1:869777..869782,+`` -> the TSS (start for +, end for -). Returns None for
    ids that don't parse (e.g. spike-ins / unmapped)."""
    m = re.match(r"(chr[\w]+):(\d+)\.\.(\d+),([+-])", row_id)
    if m:  # CAGE peak -> strand-aware TSS
        chrom, a, b, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        return chrom, (a if strand == "+" else b)
    m = re.match(r"(chr[\w]+):(\d+)-(\d+)", row_id)
    if m:  # enhancer -> midpoint
        return m.group(1), (int(m.group(2)) + int(m.group(3))) // 2
    return None


def load_gencode_tss(path: str) -> list[tuple[str, int, str, str]]:
    """Gene TSS list from a GENCODE GFF3: (chrom, tss, gene_id, gene_name)."""
    out = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            tss = start if strand == "+" else end
            gid = re.search(r"gene_id[=\s]+\"?([^\";]+)", attrs)
            gname = re.search(r"gene_name[=\s]+\"?([^\";]+)", attrs)
            if gid:
                out.append((chrom, tss, gid.group(1), gname.group(1) if gname else ""))
    return out


def assign_peaks_to_genes(peak_ids, genes, tss_window):
    """peak_index -> (gene_id, gene_name) for peaks within ``tss_window`` of a gene TSS.

    Nearest TSS wins. Peaks not near any annotated TSS are omitted (they aren't the
    promoter of a gene we can name, so they can't seed an enhancer->gene edge)."""
    by_chrom: dict[str, tuple[np.ndarray, list]] = {}
    for chrom in {g[0] for g in genes}:
        gs = sorted((g for g in genes if g[0] == chrom), key=lambda g: g[1])
        by_chrom[chrom] = (np.array([g[1] for g in gs]), gs)
    out: dict[int, tuple[str, str]] = {}
    for i, pid in enumerate(peak_ids):
        pt = parse_point(pid)
        if pt is None or pt[0] not in by_chrom:
            continue
        tss_arr, gs = by_chrom[pt[0]]
        j = int(np.searchsorted(tss_arr, pt[1]))
        best = None
        for k in (j - 1, j):  # nearest TSS is one of the two bracketing entries
            if 0 <= k < len(gs):
                d = abs(int(tss_arr[k]) - pt[1])
                if d <= tss_window and (best is None or d < best[0]):
                    best = (d, gs[k][2], gs[k][3])
        if best:
            out[i] = (best[1], best[2])
    return out


def bh_fdr(pvals: np.ndarray, m: int) -> np.ndarray:
    """Benjamini-Hochberg q-values for a subset of size len(pvals) drawn from ``m``
    total tests. Using the true ``m`` (not len(pvals)) keeps FDR correct even though
    we only stored the r>R_FLOOR pairs."""
    order = np.argsort(pvals)
    ranked = pvals[order] * m / (np.arange(len(pvals)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonic step-up
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def generate(enh_matrix, peak_matrix, gencode, facets, out, *, window, tss_window,
             r_floor, fdr, id_regex):
    enh_ids, enh, enh_libs = read_matrix(enh_matrix, id_regex)
    peak_ids, peak, peak_libs = read_matrix(peak_matrix, id_regex)

    if facets:
        fmap = dict(l.split("\t")[:2] for l in Path(facets).read_text().splitlines() if "\t" in l)
        enh, ef = facet_average(enh, enh_libs, fmap)
        peak, pf = facet_average(peak, peak_libs, fmap)
        # correlation requires the SAME facets in the SAME order in both matrices
        common = [f for f in ef if f in set(pf)]
        enh = enh[:, [ef.index(f) for f in common]]
        peak = peak[:, [pf.index(f) for f in common]]
    else:  # fall back to shared raw libraries (over-powered; warn)
        common_libs = [l for l in enh_libs if l in set(peak_libs)]
        print(f"WARNING: no --facets; correlating over {len(common_libs)} raw libraries "
              "-> anticonservative FDR. Provide a sample->facet map.", file=sys.stderr)
        enh = enh[:, [enh_libs.index(l) for l in common_libs]]
        peak = peak[:, [peak_libs.index(l) for l in common_libs]]

    n = enh.shape[1]
    if n < 4:
        raise SystemExit(f"only {n} shared facets/libraries -- too few to correlate")

    genes = load_gencode_tss(gencode)
    peak_gene = assign_peaks_to_genes(peak_ids, genes, tss_window)  # keep only promoter peaks
    prom_idx = np.array(sorted(peak_gene))                          # peak rows that map to a gene
    if len(prom_idx) == 0:
        raise SystemExit("no CAGE peaks fell within a gene TSS window -- check inputs/build")
    peak, peak_ids = peak[prom_idx], [peak_ids[i] for i in prom_idx]
    peak_gene = {j: peak_gene[i] for j, i in enumerate(prom_idx)}

    # z-standardize rows over facets (ddof=0) so a dot product / n is Pearson r.
    def zrows(a):
        mu = a.mean(axis=1, keepdims=True)
        sd = a.std(axis=1, ddof=0, keepdims=True)
        ok = (sd.ravel() > 0)
        z = np.zeros_like(a)
        z[ok] = (a[ok] - mu[ok]) / sd[ok]
        return z, ok

    ze, enh_ok = zrows(enh)
    zp, peak_ok = zrows(peak)

    # per-chrom sorted promoter positions for the 500 kb window lookup
    ppt = [parse_point(p) for p in peak_ids]
    chrom_of = np.array([p[0] if p else "" for p in ppt])
    pos_of = np.array([p[1] if p else -1 for p in ppt])
    chrom_idx: dict[str, np.ndarray] = {}
    for c in set(chrom_of):
        sel = np.where((chrom_of == c) & peak_ok)[0]
        chrom_idx[c] = sel[np.argsort(pos_of[sel])]

    m = 0                    # true count of tested pairs (for BH)
    rows = []                # stored candidates: (enh_i, peak_j, r)
    for ei, eid in enumerate(enh_ids):
        if not enh_ok[ei]:
            continue
        ept = parse_point(eid)
        if ept is None or ept[0] not in chrom_idx:
            continue
        chrom, mid = ept
        cand = chrom_idx[chrom]
        cpos = pos_of[cand]
        lo, hi = np.searchsorted(cpos, [mid - window, mid + window])
        sel = cand[lo:hi]
        if len(sel) == 0:
            continue
        r = (zp[sel] @ ze[ei]) / n     # vector of Pearson r's, enhancer vs each promoter
        m += len(sel)
        take = np.where(r > r_floor)[0]
        for t in take:
            rows.append((ei, int(sel[t]), float(r[t])))

    if not rows:
        raise SystemExit("no positively-correlated pairs above R_FLOOR -- nothing to emit")

    r = np.array([x[2] for x in rows])
    t = r * np.sqrt((n - 2) / np.clip(1 - r * r, 1e-12, None))
    p = 2 * stats.t.sf(np.abs(t), df=n - 2)
    q = bh_fdr(p, m)

    # collapse to one edge per (enhancer, gene), keeping the strongest correlation
    best: dict[tuple[str, str], tuple[float, float, str]] = {}
    for (ei, pj, rv), qv in zip(rows, q):
        if qv >= fdr or rv <= 0:
            continue
        gid, gname = peak_gene[pj]
        key = (enh_ids[ei], gid)
        if key not in best or rv > best[key][0]:
            best[key] = (rv, float(qv), gname)

    with open(out, "w") as fh:
        fh.write("enhancer_id\tgene_id\tgene_name\trel\tR\tFDR\n")
        for (enh_id, gid), (rv, qv, gname) in sorted(best.items()):
            fh.write(f"{enh_id}\t{gid}\t{gname}\tcoactive\t{rv:.4f}\t{qv:.3g}\n")
    print(f"wrote {len(best)} enhancer->gene coactive edges to {out} "
          f"({m} pairs tested, {len(rows)} above R_FLOOR)", file=sys.stderr)


def selftest():
    """A planted enhancer->gene correlation survives; an unrelated one doesn't."""
    n = 12
    rng = np.arange(n, dtype="float32")
    enh = np.vstack([rng, rng[::-1]])                      # enh0 ramps up, enh1 ramps down
    peak = np.vstack([rng + 0.01, np.ones(n)])             # peak0 tracks enh0; peak1 is flat
    ze, _ = (lambda a: ((a - a.mean(1, keepdims=True)) /
             np.clip(a.std(1, ddof=0, keepdims=True), 1e-9, None), None))(enh)
    zp, _ = (lambda a: ((a - a.mean(1, keepdims=True)) /
             np.clip(a.std(1, ddof=0, keepdims=True), 1e-9, None), None))(peak)
    r = (zp @ ze.T) / n
    assert r[0, 0] > 0.99, r[0, 0]        # enh0 vs peak0: near-perfect
    assert abs(r[1, 0]) < 0.99            # enh1 vs peak0: anti/decorrelated
    q = bh_fdr(np.array([1e-8, 0.5]), m=100)
    assert q[0] < 0.05 < q[1], q          # BH with true m keeps the strong hit
    print("selftest OK")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--enhancer-matrix")
    ap.add_argument("--peak-matrix")
    ap.add_argument("--gencode")
    ap.add_argument("--facets", default="")
    ap.add_argument("--out", default="fantom_enhancer_gene.coactive.tsv")
    ap.add_argument("--window", type=int, default=WINDOW)
    ap.add_argument("--tss-window", type=int, default=TSS_WINDOW)
    ap.add_argument("--r-floor", type=float, default=R_FLOOR)
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--id-regex", default=ID_RE)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    for req in ("enhancer_matrix", "peak_matrix", "gencode"):
        if not getattr(a, req):
            ap.error(f"--{req.replace('_', '-')} is required")
    generate(a.enhancer_matrix, a.peak_matrix, a.gencode, a.facets, a.out,
             window=a.window, tss_window=a.tss_window, r_floor=a.r_floor,
             fdr=a.fdr, id_regex=a.id_regex)


if __name__ == "__main__":
    main()
