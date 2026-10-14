#!/usr/bin/env python3
"""
process_crispri.py — 10x CRISPRi screen  ->  processed .h5ad for the site figure

Takes a Cell Ranger `filtered_feature_bc_matrix.h5` plus a geomux guide-assignment
table and produces an AnnData with PCA, UMAP, cell-cycle phase, and a per-cell
heterogeneity score, ready for `tools/export_embedding.py`.

    python3 tools/process_crispri.py \
        --h5     .../filtered_feature_bc_matrix.h5 \
        --guides .../geomux_output_assignments_CRISPRi.csv \
        --output processed.h5ad

It also reproduces the heterogeneity analysis of Figure 3C-D: genome-wide
transcript coefficient of variation per perturbation, each compared against a
*size-matched* random sample of non-targeting control cells. Size matching is
not optional — dispersion estimates scale with the number of cells, so unequal
group sizes alone manufacture an apparent heterogeneity difference.

CV is computed on linear relative counts, never on log-transformed values: the
CV of a log is a different quantity, and taking it on log data flattens the
effect to nothing.

Only standard scanpy steps are used; every threshold is a flag so the choices
stay visible rather than buried.
"""

from __future__ import annotations

import argparse
import json
import re
import sys

# Regev-lab cell-cycle genes, as distributed with Seurat and scanpy tutorials.
S_GENES = """MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1
HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51
RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1
E2F8""".split()

G2M_GENES = """HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B
MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E
TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5
CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3
GAS2L3 CBX5 CENPA""".split()


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5", required=True, help="Cell Ranger filtered_feature_bc_matrix.h5")
    p.add_argument("--guides", required=True, help="geomux assignment CSV (needs barcode + guide)")
    p.add_argument("--output", required=True, help="path for the processed .h5ad")
    p.add_argument("--min-genes", type=int, default=200, help="drop cells below this many genes")
    p.add_argument("--min-cells", type=int, default=3, help="drop genes seen in fewer cells")
    p.add_argument("--max-mt", type=float, default=20.0, help="drop cells above this %% mitochondrial")
    p.add_argument("--n-hvg", type=int, default=2000, help="highly variable genes for PCA")
    p.add_argument("--n-pcs", type=int, default=50, help="principal components to compute")
    p.add_argument("--n-pcs-metric", type=int, default=30,
                   help="PCs used for the neighbor graph")
    p.add_argument("--min-cells-per-group", type=int, default=50,
                   help="perturbation groups smaller than this are dropped (a centroid "
                        "from a handful of cells gives a meaningless CV)")
    p.add_argument("--keep-targets", default="",
                   help="comma-separated targets to retain for DISPLAY (e.g. NT,RNF8,MIS18A). "
                        "CV statistics are still computed on every group beforehand.")
    p.add_argument("--balance", action="store_true",
                   help="downsample each retained group to the smallest, so the figure is "
                        "not showing a group-size difference as a density difference")
    p.add_argument("--control", default="NT", help="target name of the non-targeting control")
    p.add_argument("--cv-cutoff", type=float, default=0.1,
                   help="gene filter: median relative counts must reach this in both groups")
    p.add_argument("--cv-reps", type=int, default=20,
                   help="random size-matched control draws to average over")
    p.add_argument("--draw-reps", type=int, default=200,
                   help="random equal-size draws used to report how often each group "
                        "is tighter than the control")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    import numpy as np
    import pandas as pd
    import scanpy as sc

    sc.settings.verbosity = 1
    np.random.seed(args.seed)

    def say(msg):
        print(f"[process] {msg}", file=sys.stderr, flush=True)

    # --- load ----------------------------------------------------------------
    say(f"reading {args.h5}")
    adata = sc.read_10x_h5(args.h5, gex_only=False)
    adata.var_names_make_unique()
    say(f"  {adata.n_obs:,} cells x {adata.n_vars:,} features")

    if "feature_types" in adata.var:
        keep = adata.var["feature_types"] == "Gene Expression"
        say(f"  keeping {int(keep.sum()):,} Gene Expression features "
            f"(dropping {int((~keep).sum())} guide-capture)")
        adata = adata[:, keep].copy()

    # --- guide assignments ---------------------------------------------------
    say(f"reading {args.guides}")
    g = pd.read_csv(args.guides)
    # The geomux output repeats the barcode column; take the first occurrence.
    g = g.loc[:, ~g.columns.duplicated()]
    if "barcode" not in g.columns or "guide" not in g.columns:
        sys.exit(f"guide table needs 'barcode' and 'guide' columns; found {list(g.columns)}")
    g = g[["barcode", "guide"]].drop_duplicates("barcode").set_index("barcode")

    overlap = adata.obs_names.intersection(g.index)
    say(f"  {len(g):,} assignments, {len(overlap):,} match the matrix "
        f"({100*len(overlap)/len(g):.1f}%)")
    if len(overlap) < 100:
        sys.exit("too few barcodes matched — check that the guide table and matrix "
                 "come from the same run")

    adata = adata[overlap].copy()
    adata.obs["guide"] = g.loc[adata.obs_names, "guide"].astype(str).values

    # Collapse "MIS18A-1i" / "MIS18A-2i" to the target gene. Guides against the
    # same gene are replicates of one perturbation, not separate conditions.
    def target_of(guide: str) -> str:
        m = re.match(r"^(.+?)-\d+[ia]$", guide)
        return m.group(1) if m else guide

    _t = [target_of(x) for x in adata.obs["guide"]]
    # Control first, then alphabetical: the figure colors by level index,
    # so this ordering is the legend's ordering.
    adata.obs["target"] = pd.Categorical(
        _t, categories=[args.control] + sorted(set(_t) - {args.control}), ordered=True)
    # A 4-level view for the figure legend: the control, the two hits, everything else.
    _pert = ["Non-targeting" if t == args.control else t if t in ("RNF8", "MIS18A")
             else "Other chromatin organizer" for t in adata.obs["target"]]
    # Ordered, not alphabetical: control first, then the two hits. The figure
    # assigns colors by level index, so this order is the legend's order.
    adata.obs["perturbation"] = pd.Categorical(
        _pert, categories=["Non-targeting", "RNF8", "MIS18A",
                           "Other chromatin organizer"], ordered=True)
    say(f"  {adata.obs['target'].nunique()} target genes, "
        f"{adata.obs['guide'].nunique()} guides")

    # --- QC ------------------------------------------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mt].copy()
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    say(f"QC: {n0:,} -> {adata.n_obs:,} cells "
        f"({n0 - adata.n_obs} dropped), {adata.n_vars:,} genes retained")

    # Drop perturbation groups too small for a meaningful CV estimate.
    counts = adata.obs["target"].value_counts()
    small = counts[counts < args.min_cells_per_group]
    if len(small):
        say(f"dropping {len(small)} group(s) under {args.min_cells_per_group} cells: "
            f"{', '.join(small.index)}")
        adata = adata[adata.obs["target"].isin(counts[counts >= args.min_cells_per_group].index)].copy()

    # --- normalize -----------------------------------------------------------
    adata.layers["counts"] = adata.X.copy()
    full_counts = adata.X.copy()   # all genes, pre-HVG — the CV analysis needs these
    control_name = args.control
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # Cell cycle is scored on all genes, before HVG subsetting.
    s = [x for x in S_GENES if x in adata.var_names]
    g2m = [x for x in G2M_GENES if x in adata.var_names]
    say(f"cell cycle: {len(s)}/{len(S_GENES)} S genes, {len(g2m)}/{len(G2M_GENES)} G2M genes found")
    sc.tl.score_genes_cell_cycle(adata, s_genes=s, g2m_genes=g2m)

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    # --- embed ---------------------------------------------------------------
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs_metric, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    say("PCA + UMAP done")

    # --- heterogeneity: genome-wide CV vs size-matched controls --------------
    # The published metric (Fig. 3C-D). Each perturbation is compared against a
    # random NT sample of the SAME size, because CV estimates scale with n.
    from scipy import sparse

    counts = full_counts
    tot = np.asarray(counts.sum(axis=1)).ravel()
    RC = (sparse.diags(1e4 / np.maximum(tot, 1)) @ counts).tocsr()

    targets_arr = adata.obs["target"].values
    nt_idx = np.where(targets_arr == control_name)[0]
    rng = np.random.default_rng(args.seed)

    def cv_pair(sel_t, sel_c, cutoff):
        """Median CV of each group over genes detectable in both."""
        T = RC[sel_t].toarray()
        C = RC[sel_c].toarray()
        det = (np.median(T, axis=0) >= cutoff) & (np.median(C, axis=0) >= cutoff)
        if det.sum() < 50:
            return None
        def cv(a):
            m = a[:, det].mean(0)
            sd = a[:, det].std(0)
            return np.where(m > 0, sd / m, np.nan)
        ct, cc = cv(T), cv(C)
        ok = np.isfinite(ct) & np.isfinite(cc)
        return dict(
            n_genes=int(det.sum()),
            cv=float(np.nanmedian(ct)),
            cv_ratio=float(np.nanmedian(ct) / np.nanmedian(cc)),
            pct_down=float(100 * np.mean(ct[ok] < cc[ok])),
        )

    say(f"CV vs size-matched {control_name} controls "
        f"({args.cv_reps} control draws, gene cutoff {args.cv_cutoff}):")
    cv_stats = {}
    for t in sorted(set(targets_arr)):
        idx = np.where(targets_arr == t)[0]
        n = len(idx)
        if t == control_name:
            # Split NT against itself: the method's own null. Should sit at 1.0.
            reps = []
            for _ in range(args.cv_reps):
                perm = rng.permutation(nt_idx)
                half = min(n // 2, len(perm) // 2)
                r = cv_pair(perm[:half], perm[half:2 * half], args.cv_cutoff)
                if r: reps.append(r)
        else:
            pool = np.setdiff1d(nt_idx, idx)
            if len(pool) < n:
                say(f"    {t}: {n} cells exceeds the {len(pool)}-cell control pool; skipped")
                continue
            reps = []
            for _ in range(args.cv_reps):
                r = cv_pair(idx, rng.choice(pool, n, replace=False), args.cv_cutoff)
                if r: reps.append(r)
        if not reps:
            continue
        cv_stats[t] = dict(
            n_cells=int(n),
            n_genes=int(np.median([r["n_genes"] for r in reps])),
            cv_ratio=round(float(np.mean([r["cv_ratio"] for r in reps])), 4),
            cv_ratio_sd=round(float(np.std([r["cv_ratio"] for r in reps])), 4),
            pct_down=round(float(np.mean([r["pct_down"] for r in reps])), 1),
        )

    for t, v in sorted(cv_stats.items(), key=lambda kv: kv[1]["cv_ratio"]):
        mark = "  <- control (null)" if t == control_name else (
            "  **" if t in ("RNF8", "MIS18A") else "")
        say(f"    {t:<10} n={v['n_cells']:>5}  CV ratio {v['cv_ratio']:.3f}"
            f" +/-{v['cv_ratio_sd']:.3f}   {v['pct_down']:>5.1f}% genes down{mark}")

    adata.uns["cv_stats"] = cv_stats

    # --- the paper's two specificity controls --------------------------------
    # 1. Cell cycle: the effect should not be a cycling artifact.
    ph = pd.crosstab(adata.obs["target"], adata.obs["phase"], normalize="index") * 100
    say("cell-cycle composition (%), control vs the two hits:")
    for t in [control_name, "RNF8", "MIS18A"]:
        if t in ph.index:
            say(f"    {t:<10} " + "  ".join(f"{c} {ph.loc[t, c]:5.1f}" for c in ph.columns))

    # 2. Mean expression: a CV shift must not be a mean shift in disguise.
    #    Compared per gene — the global mean of relative counts is identical by
    #    construction, so quoting that would prove nothing.
    say("per-gene mean expression vs control (median ratio; 1.0 = no mean shift):")
    nt_mean = np.asarray(RC[nt_idx].mean(axis=0)).ravel()
    for t in [control_name, "RNF8", "MIS18A"]:
        idx = np.where(targets_arr == t)[0]
        if not len(idx):
            continue
        tm = np.asarray(RC[idx].mean(axis=0)).ravel()
        det = (nt_mean > 0.05) & (tm > 0.05)
        ratio = np.median(tm[det] / nt_mean[det])
        rho = float(np.corrcoef(np.log1p(tm[det]), np.log1p(nt_mean[det]))[0, 1])
        say(f"    {t:<10} median ratio {ratio:.3f}   log-mean r = {rho:.4f}  "
            f"({int(det.sum())} genes)")

    # Snapshot the full cohort before subsetting — the draw-frequency statistic
    # below needs groups at their true sizes.
    pre_pca = adata.obsm["X_pca"].copy()
    pre_umap = adata.obsm["X_umap"].copy()
    pre_target = np.asarray(adata.obs["target"].values)

    # --- display subset ------------------------------------------------------
    # Done last, after every statistic above, so the numbers come from all the
    # cells that exist while the figure shows a balanced view of them.
    keep = [x.strip() for x in args.keep_targets.split(",") if x.strip()]
    if keep:
        missing = [k for k in keep if k not in set(adata.obs["target"])]
        if missing:
            sys.exit(f"--keep-targets names absent group(s): {', '.join(missing)}")
        adata = adata[adata.obs["target"].isin(keep)].copy()
        say(f"display subset: {', '.join(keep)} ({adata.n_obs:,} cells)")

    if args.balance:
        counts = adata.obs["target"].value_counts()
        counts = counts[counts > 0]
        n_each = int(counts.min())
        rng2 = np.random.default_rng(args.seed)
        take = np.concatenate([
            rng2.choice(np.where(adata.obs["target"].values == g)[0], n_each, replace=False)
            if counts[g] > n_each else np.where(adata.obs["target"].values == g)[0]
            for g in counts.index])
        adata = adata[np.sort(take)].copy()
        say(f"balanced to {n_each:,} cells per group (seed {args.seed}); "
            f"{adata.n_obs:,} shown")
        adata.uns["display"] = {
            "balanced": True, "cellsPerGroup": n_each, "seed": int(args.seed),
            "groups": list(counts.index),
        }

        # How often does a random equal-size draw put each group below the
        # control? Measured on the FULL groups captured before balancing —
        # drawing n cells from a group that has exactly n is not sampling, and
        # would report a meaningless 100%.
        PCf = pre_pca[:, : args.n_pcs_metric]
        freq = {}
        for space, X in (("pc", PCf), ("umap", pre_umap)):
            for g in counts.index:
                if g == control_name:
                    continue
                pool_g = np.where(pre_target == g)[0]
                pool_c = np.where(pre_target == control_name)[0]
                wins = 0
                for sdraw in range(args.draw_reps):
                    r3 = np.random.default_rng(sdraw)
                    a = r3.choice(pool_g, min(n_each, len(pool_g)), replace=False)
                    b = r3.choice(pool_c, min(n_each, len(pool_c)), replace=False)
                    sp = lambda sel: np.median(
                        np.linalg.norm(X[sel] - X[sel].mean(0), axis=1))
                    wins += sp(a) < sp(b)
                freq[f"{g}_{space}"] = round(100 * wins / args.draw_reps, 1)
        adata.uns["display"]["tighterThanControlPct"] = freq
        adata.uns["display"]["drawReps"] = int(args.draw_reps)
        say(f"tighter than control across {args.draw_reps} random equal-size draws "
            f"of the full groups:")
        for k, v in sorted(freq.items()):
            say(f"    {k:<16} {v:5.1f}%")

    # Filtering upstream resets categorical ordering, so set it here — last —
    # where it will survive to the export.
    if hasattr(adata.obs["target"], "cat"):
        cats = list(adata.obs["target"].cat.categories)
        # Follow the order given on the command line; otherwise control first.
        wanted = [k for k in keep if k in cats] if keep else []
        order = wanted or ([control_name] + [c for c in cats if c != control_name])
        order += [c for c in cats if c not in order]
        adata.obs["target"] = adata.obs["target"].cat.reorder_categories(order)

    adata.write_h5ad(args.output)
    say(f"wrote {args.output}  ({adata.n_obs:,} cells, {adata.n_vars:,} HVGs)")

    summary = {
        "cells": int(adata.n_obs),
        "targets": int(adata.obs["target"].nunique()),
        "guides": int(adata.obs["guide"].nunique()),
        "cv_stats": cv_stats,
    }
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
