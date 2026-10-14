# ssobt.github.io

Personal site for Simon Sobti — computational biologist, PhD UCSF.

Hand-written HTML, CSS, and JavaScript. **No build step, no framework, no
dependencies, no webfonts.** Clone it, open it, edit it; it will still work in
five years.

## Run it locally

```bash
python3 -m http.server 8000
# open http://localhost:8000
```

A server is required (not `file://`) — the page uses ES modules and `fetch`,
both of which are blocked under the `file:` origin.

## Layout

```
index.html                  single-page main site
projects/                   per-project deep dives
assets/css/site.css         design tokens + all styling
assets/js/scatter.js        canvas embedding renderer (no charting library)
assets/js/cohort.js         patient-cohort figure: PC space + distance strip
tools/process_crispri.py    10x CRISPRi screen -> processed .h5ad
assets/js/explorer.js       controls, and the boot for both figures
assets/js/site.js           page chrome: scrollspy, reveal, copy-email
assets/data/embedding.json  the cells the single-cell viewer draws
assets/data/cohort.json     the patients the cohort figure draws
assets/data/SCHEMA.md       format contract for both
tools/                      data export, validation, image generation
```

## The data

Both figures render **real data**, regenerated from source with the tools below.

| Figure | Data | Built by |
|---|---|---|
| §02 step 1 | 1,100 TCGA-BRCA tumors (RNA-Seq V2 RSEM, cBioPortal) | `tools/export_cohort.R` |
| §02 step 3 | 7,474 MDA-MB-231 CRISPRi Perturb-seq cells | `tools/process_crispri.py` → `tools/export_embedding.py` |

### Patient cohort

```bash
Rscript tools/export_cohort.R \
  --input      /path/to/data_RNA_Seq_v2_expression_median.txt \
  --genes      RNF8,MIS18A \
  --controls   ADORA2B,ARHGEF5 \
  --pcs        50 \
  --null-genes 200 \
  --cohort     "TCGA-BRCA"
```

Reads either an `.rds` matrix or a cBioPortal `.txt` directly. Beyond the
quartile split it builds an **empirical null** from randomly chosen expressed
genes, because at this cohort size ~84% of arbitrary genes reach `P < 0.05` on
the spread metric — the z-score against that null is the meaningful number, and
the figure reports it instead of the p-value.

### Single-cell screen

```bash
# raw 10x matrix + guide assignments -> processed .h5ad
python3 tools/process_crispri.py \
  --h5     /path/to/filtered_feature_bc_matrix.h5 \
  --guides /path/to/geomux_output_assignments_CRISPRi.csv \
  --output crispri.h5ad

# .h5ad -> embedding.json
python3 tools/export_embedding.py \
  --input crispri.h5ad --output assets/data/embedding.json \
  --reductions X_pca,X_umap --cloud X_pca \
  --labels perturbation,target,phase \
  --group-stats cv_stats --group-stats-for target

python3 tools/validate_embedding.py assets/data/embedding.json
```

`process_crispri.py` also reproduces the paper's Figure 3C–D: genome-wide
transcript CV per perturbation against a **size-matched** non-targeting control
sample, on linear relative counts. Both details matter — dispersion scales with
cell number, and the CV of a log is a different quantity.

Needs scanpy; the environment is not one of the repo's dependencies:

```bash
mamba create -n sc_export -c conda-forge python=3.11 scanpy anndata h5py \
  numpy pandas scipy scikit-learn umap-learn python-igraph leidenalg
```

### Synthetic fallbacks

`tools/make_placeholder.py` and `tools/make_cohort_placeholder.py` generate
synthetic stand-ins with the same schema, for working on the figures without the
source data to hand. Their output labels itself as simulated in the caption.

The exporters emit **only** coordinates and the columns you name — never counts,
barcodes, or patient identifiers.

> TCGA gene expression from the GDC is **open access** (the controlled tier is
> raw sequence and germline variants), so publishing PCA coordinates derived
> from it is fine. Just never emit patient barcodes — the exporters write
> anonymous indices.

## Regenerating images

```bash
python3 tools/make_placeholder.py          # synthetic embedding.json
python3 tools/make_cohort_placeholder.py   # synthetic cohort.json
python3 tools/make_og_image.py             # assets/img/og.svg
# then rasterize og.svg -> og.png with headless Chrome (see that file's docstring)
```

## Deploying

The repo is a GitHub user site: pushing to `main` publishes the root of the repo
at `https://ssobt.github.io`. `.nojekyll` keeps Jekyll out of the way.

## Notes for future edits

- Header and footer markup is duplicated across the three HTML files on purpose.
  With no build step, duplication is more robust than a runtime include, and it
  keeps the pages working with JavaScript disabled.
- All color, spacing, and type decisions live as custom properties at the top of
  `site.css`. Change them there, not inline.
- The categorical palette is Okabe–Ito, chosen to stay distinguishable under the
  common forms of color vision deficiency. If you add series, extend
  `CATEGORICAL` in `scatter.js` rather than hard-coding colors.
