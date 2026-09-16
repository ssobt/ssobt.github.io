#!/usr/bin/env python3
"""
make_cohort_placeholder.py — generate a SYNTHETIC cohort.json

Stands in for the TCGA-BRCA analysis in Figure 1 of Woo & Sobti et al. (2026)
until real coordinates are exported with tools/export_cohort.R. The output is
labelled synthetic in `meta`, which the page renders into the figure legend.

The simulation reproduces the *structure* of the published result, not its
values:

  - 1,200 patients live in a 50-dimensional latent space; dimensions 1-2 are
    shown as PC1/PC2, exactly as the paper visualises PC space while computing
    its spread metric across the top 50 PCs.
  - For a *candidate* regulator, expression is coupled to a patient's radial
    position, so top-quartile expressers sit further from their group centroid.
    That is the published finding: high CO expression tracks with greater
    inter-tumour transcriptomic spread.
  - For a *control* gene, expression is independent of position, so the two
    quartiles have the same spread. Toggling to a control is what shows the
    effect is specific rather than an artifact of quartile-splitting anything.

Usage:
    python3 tools/make_cohort_placeholder.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random

SEED = 20260916
N = 1200
DIMS = 50

# Per-dimension standard deviations, falling off like a real PCA scree plot.
def dim_sd(i: int) -> float:
    return 12.0 * math.exp(-i / 14.0) + 1.1


# Published values from Woo & Sobti et al., bioRxiv 2026 (doi:10.64898/2026.04.18.719392).
# `centroidP` is Figure 1D (two-sided Wilcoxon signed-rank on distance from
# quartile centroid); `hr` / `survP` are Figure 6C (METABRIC, log-rank).
GENES = [
    {"key": "RNF8",     "role": "candidate", "centroidP": "3e-12",   "hr": 1.24, "survP": "0.01"},
    {"key": "MIS18A",   "role": "candidate", "centroidP": "5.8e-11", "hr": 1.75, "survP": "<0.001"},
    {"key": "ADORA2B",  "role": "control",   "centroidP": "0.17",    "hr": 1.09, "survP": "0.34"},
    {"key": "ARHGEF5",  "role": "control",   "centroidP": "0.5",     "hr": 0.72, "survP": "<0.001"},
]

Q25, MID, Q75 = 0, 1, 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", default="assets/data/cohort.json")
    ap.add_argument("--n", type=int, default=N)
    args = ap.parse_args()

    rng = random.Random(SEED)
    n = args.n

    # --- latent space -------------------------------------------------------
    # Radial scale varies per patient so the cloud has a genuine spread
    # gradient for candidate expression to couple to.
    scale = [math.exp(rng.gauss(0, 0.42)) for _ in range(n)]
    pts = [[rng.gauss(0, dim_sd(d)) * scale[i] for d in range(DIMS)] for i in range(n)]

    # Radial position in the full latent space, standardised.
    radius = [math.sqrt(sum(v * v for v in p)) for p in pts]
    r_mean = sum(radius) / n
    r_sd = math.sqrt(sum((r - r_mean) ** 2 for r in radius) / n)
    r_z = [(r - r_mean) / r_sd for r in radius]

    def centroid(idx: list[int]) -> list[float]:
        return [sum(pts[i][d] for i in idx) / len(idx) for d in range(DIMS)]

    def dist_to(p: list[float], c: list[float]) -> float:
        return math.sqrt(sum((p[d] - c[d]) ** 2 for d in range(DIMS)))

    genes = {}
    for g in GENES:
        # Candidates: expression tracks radial position (the published effect).
        # Controls: expression is independent of where the patient sits.
        if g["role"] == "candidate":
            latent = [0.62 * r_z[i] + rng.gauss(0, 0.78) for i in range(n)]
        else:
            latent = [rng.gauss(0, 1) for _ in range(n)]

        # FPKM-like values, for the tooltip only.
        expr = [round(math.exp(2.4 + 0.55 * z), 2) for z in latent]

        order = sorted(range(n), key=lambda i: latent[i])
        cut = n // 4
        quart = [MID] * n
        for i in order[:cut]:
            quart[i] = Q25
        for i in order[-cut:]:
            quart[i] = Q75

        q25_idx = [i for i in range(n) if quart[i] == Q25]
        q75_idx = [i for i in range(n) if quart[i] == Q75]

        # Distance to *that group's own* centroid, across all 50 dimensions —
        # the metric defined in the paper's methods.
        c25, c75 = centroid(q25_idx), centroid(q75_idx)
        dist = [0.0] * n
        for i in q25_idx:
            dist[i] = dist_to(pts[i], c25)
        for i in q75_idx:
            dist[i] = dist_to(pts[i], c75)
        # Middle patients are not part of either comparison; measured against
        # the cohort centroid purely so the tooltip has something to show.
        mid_idx = [i for i in range(n) if quart[i] == MID]
        cmid = centroid(mid_idx)
        for i in mid_idx:
            dist[i] = dist_to(pts[i], cmid)

        genes[g["key"]] = {
            "label": g["key"],
            "role": g["role"],
            "expression": expr,
            "quartile": quart,
            "distance": [round(d, 2) for d in dist],
            "published": {
                "centroidP": g["centroidP"],
                "hr": g["hr"],
                "survP": g["survP"],
            },
        }

    coords = []
    for p in pts:
        coords += [round(p[0], 3), round(p[1], 3)]

    payload = {
        "meta": {
            "n": n,
            "cohort": "TCGA-BRCA (simulated)",
            "synthetic": True,
            "seed": SEED,
            "dims": DIMS,
            "caption": (
                f"{n:,} simulated breast tumours, standing in for the TCGA-BRCA cohort. "
                "Patients are split into bottom (q25) and top (q75) quartiles of the "
                "selected gene's expression, and each patient's distance to its own "
                "group's centroid is measured across the top 50 principal components. "
                "For a candidate regulator the top quartile is measurably more spread "
                "out; for a negative control the two quartiles look the same. This is "
                "simulated data reproducing the structure of the published result, not "
                "the result itself — the real statistics are cited beside each gene."
            ),
            "method": (
                "Patients projected into PCA space; spread measured as distance from "
                "each patient to its quartile group's centroid across the top 50 PCs."
            ),
            "published": {
                "citation": "Woo BJ*, Sobti S*, et al. bioRxiv 2026, Figure 1C–D and 6C",
                "doi": "10.64898/2026.04.18.719392",
            },
            "exportedBy": "tools/make_cohort_placeholder.py",
        },
        "coords": {"pc": coords},
        "geneOrder": [g["key"] for g in GENES],
        "genes": genes,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    kb = round(os.path.getsize(args.output) / 1024)
    print(f"Wrote {args.output} ({kb} KB, {n} synthetic patients).")
    for g in GENES:
        gg = genes[g["key"]]
        d25 = sorted(gg["distance"][i] for i in range(n) if gg["quartile"][i] == Q25)
        d75 = sorted(gg["distance"][i] for i in range(n) if gg["quartile"][i] == Q75)
        m = lambda a: a[len(a) // 2]
        print(f"  {g['key']:<9} {g['role']:<10} median q25={m(d25):6.2f}  "
              f"q75={m(d75):6.2f}  ratio={m(d75)/m(d25):.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
