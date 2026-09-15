#!/usr/bin/env python3
"""
export_embedding.py — AnnData (scanpy)  ->  assets/data/embedding.json

Emits ONLY obsm coordinates and the obs columns you name. The counts matrix,
barcodes, and every unnamed obs column stay behind.

Usage
    python3 tools/export_embedding.py \
        --input      adata.h5ad \
        --output     assets/data/embedding.json \
        --reductions X_pca,X_umap \
        --labels     guide,phase \
        --scores     het_cv \
        --n          6000 \
        --source     "MDA-MB-231 CRISPRi Perturb-seq" \
        --caption    "6,000 cells from ..."

Only --input is required.

Requires anndata (and numpy/pandas, which come with it):
    pip install anndata
"""

from __future__ import annotations

import argparse
import json
import os
import sys

# Common obsm keys mapped to the short names the site uses in its UI.
PRETTY = {"X_pca": "pca", "X_umap": "umap", "X_tsne": "tsne", "X_diffmap": "diffmap"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Export an AnnData embedding to the site's embedding.json format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, help="path to a .h5ad file")
    p.add_argument("--output", default="assets/data/embedding.json")
    p.add_argument("--reductions", default="X_pca,X_umap",
                   help="comma-separated obsm keys; order matters, the last is the default view")
    p.add_argument("--cloud", default="",
                   help="obsm key for the hero's 3D intro cloud (dims 1-3); "
                        "defaults to the first --reductions entry")
    p.add_argument("--labels", default="", help="comma-separated categorical obs columns")
    p.add_argument("--scores", default="", help="comma-separated continuous obs columns")
    p.add_argument("--n", type=int, default=6000, help="max cells to export")
    p.add_argument("--seed", type=int, default=42, help="subsampling seed")
    p.add_argument("--source", default="", help="short description of the dataset")
    p.add_argument("--caption", default="", help="figure legend; auto-generated if omitted")
    p.add_argument("--digits", type=int, default=3, help="coordinate rounding")
    return p.parse_args(argv)


def split_csv(s: str) -> list[str]:
    return [t.strip() for t in s.split(",") if t.strip()]


def main(argv=None) -> int:
    args = parse_args(argv)

    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
    except ImportError as e:
        sys.exit(f"missing dependency: {e.name}. Try:  pip install anndata")

    print(f"Reading {args.input} ...", file=sys.stderr)
    adata = ad.read_h5ad(args.input)

    reductions = split_csv(args.reductions)
    label_cols = split_csv(args.labels)
    score_cols = split_csv(args.scores)

    missing = [r for r in reductions if r not in adata.obsm]
    if missing:
        sys.exit(f"obsm key(s) not found: {', '.join(missing)}\n"
                 f"  available: {', '.join(adata.obsm.keys())}")

    missing = [c for c in label_cols + score_cols if c not in adata.obs.columns]
    if missing:
        sys.exit(f"obs column(s) not found: {', '.join(missing)}\n"
                 f"  available: {', '.join(adata.obs.columns)}")

    # --- deterministic subsample --------------------------------------------
    n_total = adata.n_obs
    rng = np.random.default_rng(args.seed)
    if n_total > args.n:
        keep = np.sort(rng.choice(n_total, size=args.n, replace=False))
    else:
        keep = np.arange(n_total)
    n = int(keep.size)
    print(f"Exporting {n} of {n_total} cells (seed {args.seed}).", file=sys.stderr)

    # --- 3D cloud for the hero intro ----------------------------------------
    cloud_key = args.cloud or reductions[0]
    coords3d = None
    if cloud_key in adata.obsm:
        m = np.asarray(adata.obsm[cloud_key])
        if m.shape[1] >= 3:
            coords3d = np.round(m[keep, :3].astype(float), args.digits).ravel().tolist()
            print(f"3D hero cloud from '{cloud_key}' dims 1-3.", file=sys.stderr)
        else:
            print(f"NOTE: '{cloud_key}' has <3 dimensions; hero will start flat.",
                  file=sys.stderr)

    # --- coordinates ---------------------------------------------------------
    coords = {}
    for key in reductions:
        m = np.asarray(adata.obsm[key])
        if m.shape[1] < 2:
            sys.exit(f"obsm['{key}'] has fewer than 2 dimensions")
        flat = np.round(m[keep, :2].astype(float), args.digits).ravel()
        coords[PRETTY.get(key, key.lstrip("X_"))] = flat.tolist()

    short = [PRETTY.get(k, k.lstrip("X_")) for k in reductions]

    # --- categorical label sets ---------------------------------------------
    label_sets = {}
    for col in label_cols:
        s = adata.obs[col].iloc[keep]
        cat = s.astype("category").cat.remove_unused_categories()
        levels = [str(x) for x in cat.cat.categories]
        if len(levels) > 24:
            print(f"WARNING: '{col}' has {len(levels)} levels; the legend will be "
                  f"unwieldy. Consider collapsing it.", file=sys.stderr)
        label_sets[col] = {
            "label": col,
            "levels": levels,
            "values": [int(v) for v in cat.cat.codes],  # already zero-indexed
        }

    # --- continuous scores ---------------------------------------------------
    scores = {}
    for col in score_cols:
        v = pd.to_numeric(adata.obs[col].iloc[keep], errors="coerce").to_numpy(dtype=float)
        if np.isnan(v).any():
            print(f"WARNING: '{col}' contains NaN; replacing with the column median.",
                  file=sys.stderr)
            v = np.nan_to_num(v, nan=float(np.nanmedian(v)))
        v = np.round(v, args.digits)
        scores[col] = {
            "label": col,
            "range": [float(v.min()), float(v.max())],
            "values": v.tolist(),
        }

    # --- assemble ------------------------------------------------------------
    caption = args.caption or (
        f"{n:,} cells. {short[-1].upper()} embedding; colour by "
        f"{label_cols[0] if label_cols else 'cluster'}."
    )

    payload = {
        "meta": {
            "n": n,
            "source": args.source,
            "caption": caption,
            "heroCaption": f"{n:,} cells · {short[-1].upper()}",
            "embeddingLabels": {s: s.upper() for s in short},
            "exportedBy": "tools/export_embedding.py",
            "seed": args.seed,
        },
        "coords": coords,
        "coords3d": coords3d,
        "labelSets": label_sets,
        "scores": scores,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    kb = round(os.path.getsize(args.output) / 1024)
    print(f"Wrote {args.output} ({kb} KB, {n} cells).", file=sys.stderr)
    if kb > 1024:
        print("NOTE: over 1 MB. Lower --n so the page still paints quickly.", file=sys.stderr)
    print(f"Validate with:  python3 tools/validate_embedding.py {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
