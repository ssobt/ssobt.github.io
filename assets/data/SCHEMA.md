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
- No patient-derived data. Use **cell-line** data (Perturb-seq, CRISPRi, scATAC).
  Individual-level TCGA expression is dbGaP-controlled; a public web page is
  not a permitted destination for it.
- Label levels are names you are willing to publish (`RNF8-Ci` is fine;
  `patient_04_pretreatment` is not).

## Validating a file

```bash
python3 tools/validate_embedding.py assets/data/embedding.json
```
