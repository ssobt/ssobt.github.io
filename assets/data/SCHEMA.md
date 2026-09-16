# `embedding.json` — format contract

The explorer in §02 of the site reads exactly one file: `assets/data/embedding.json`.
Anything conforming to this schema works with no code changes. Swap the file, reload.

Produce it with `tools/export_embedding.R` (Seurat) or `tools/export_embedding.py`
(AnnData/scanpy) rather than by hand.

---

## Shape

```jsonc
{
  "meta": {
    "n": 6000,
    "source": "MDA-MB-231 CRISPRi Perturb-seq",
    "caption": "6,000 cells …",          // rendered as the figure legend
    "heroCaption": "6,000 cells · UMAP", // short label under the hero panel
    "embeddingLabels": { "pca": "PCA", "umap": "UMAP" }
  },

  // One entry per embedding. Flat, interleaved [x0,y0,x1,y1,…].
  // Every array must have the same length: 2 × n.
  "coords": {
    "pca":  [1.234, -0.881, …],
    "umap": [8.112,  3.407, …]
  },

  // OPTIONAL. Flat [x0,y0,z0,…] of length 3 × n — normally PC1/PC2/PC3.
  // Used only by the hero's 3D intro, which rotates this cloud and then
  // collapses it onto the 2D embedding. Omit it and the hero starts flat.
  "coords3d": [1.234, -0.881, 0.402, …],

  // Categorical colourings. `values[i]` indexes into `levels`.
  "labelSets": {
    "perturbation": {
      "label":  "Perturbation",          // shown in the UI; defaults to the key
      "levels": ["NT", "RNF8-Ci", "MIS18A-Ci"],
      "values": [0, 0, 1, 2, …]          // length n, integers 0..levels.length-1
    }
  },

  // Continuous colourings, rendered on a viridis ramp with a colorbar.
  "scores": {
    "heterogeneity_cv": {
      "label":  "Transcriptional heterogeneity (CV)",
      "range":  [0.0, 1.2],              // optional; inferred from data if absent
      "values": [0.41, 0.63, …]          // length n
    }
  }
}
```

## Rules

| Rule | Why |
|---|---|
| `coords` must contain **at least one** embedding. Two (`pca`, `umap`) enables the morph animation. | The PCA→UMAP transition is the figure's most memorable moment. |
| `coords3d`, if present, must be length 3 × n. A mismatch logs a warning and disables the intro rather than throwing. | The hero is decorative; a bad 3D column should never take the page down. |
| Every array in `coords`, `labelSets.*.values`, `scores.*.values` must be length-consistent with `n`. | The loader throws a readable error on mismatch rather than rendering nonsense. |
| `labelSets` values are **integer indices**, not strings. | Roughly 6× smaller on the wire than repeated label strings. |
| Round coordinates to 3 decimals. | Beyond that is noise in an embedding, and it inflates the file. |
| Keep `n` between about 3,000 and 20,000. | Below ~3k the structure looks thin; above ~20k the JSON gets large enough to hurt first paint. |
| `meta.caption` is required in practice. | It is the figure legend. An unlabelled figure on a scientist's site is a bad look. |

## Size

At 6,000 cells with two embeddings, a 3D cloud, three label sets, and one
score: **332 KB** raw, **~95 KB** gzipped — GitHub Pages gzips automatically.
If the file exceeds ~1 MB, lower `--n` in the export script rather than
dropping precision.

## What must never go in this file

The export scripts emit **only** coordinates and the label columns you name
explicitly. They never touch the counts matrix. Before publishing, confirm:

- No cell barcodes, sample IDs, or patient identifiers — not even hashed.
- **TCGA correction:** gene expression quantification from the GDC is *open
  access*; the controlled tier is raw sequence and germline variants. Publishing
  PCA coordinates derived from TCGA expression is fine — it is what the paper
  does. Just never emit TCGA barcodes alongside them.
- Label levels are names you are willing to publish (`RNF8-Ci` is fine;
  `patient_04_pretreatment` is not).

## Validating a file

```bash
python3 tools/validate_embedding.py assets/data/embedding.json
```


---

# `cohort.json` — patient-cohort figure

The figure in §02 step 01 reads `assets/data/cohort.json`. Produce it with
`tools/export_cohort.R` (real data) or `tools/make_cohort_placeholder.py`
(synthetic stand-in).

```jsonc
{
  "meta": {
    "n": 1200,
    "cohort": "TCGA-BRCA",
    "synthetic": false,
    "dims": 50,                    // PCs used for the distance metric
    "caption": "…",                // rendered as the figure legend
    "published": { "citation": "…", "doi": "10.64898/…" }
  },

  // PC1/PC2 only, flat and interleaved: [x0,y0,x1,y1,…], length 2n.
  "coords": { "pc": [12.4, -3.1, …] },

  // Order of the gene selector, left to right.
  "geneOrder": ["RNF8", "MIS18A", "ADORA2B", "ARHGEF5"],

  "genes": {
    "RNF8": {
      "label": "RNF8",
      "role": "candidate",         // "candidate" | "control"
      "expression": [ … ],         // length n, tooltip only
      "quartile":   [ 0|1|2, … ],  // 0 = bottom, 1 = middle 50%, 2 = top
      "distance":   [ … ],         // length n, distance to own group centroid
      "published": {               // shown as a cited annotation, not computed
        "centroidP": "3e-12",
        "hr": 1.24,
        "survP": "0.01"
      }
    }
  }
}
```

## Rules

| Rule | Why |
|---|---|
| Include at least one gene with `role: "control"`. | A quartile split of *any* gene yields two groups; the control is what shows the effect is specific. The figure labels controls in the selector. |
| `distance` is measured to the patient's **own quartile group's** centroid, not the cohort centroid. | That is the published metric. Measuring to a shared centroid conflates spread with group separation. |
| `published` values are **citations, not computations**. | The figure prints them under "Published result for this gene" and attributes them. If the cloud is synthetic, the live stats and the cited stats must stay visibly distinct. |
| Never emit patient barcodes. | The exporter writes anonymous indices; keep it that way. |

At 1,200 patients and four genes the file is roughly **80 KB** raw.
