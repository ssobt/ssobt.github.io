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
assets/js/explorer.js       controls bound to the renderer
assets/js/site.js           page chrome: scrollspy, reveal, copy-email
assets/data/embedding.json  the cells the viewer draws
assets/data/SCHEMA.md       format contract for the above
tools/                      data export, validation, image generation
```

## Replacing the demo data

The viewer in §02 currently renders **synthetic placeholder data**, labelled as
such in its caption. To swap in a real export:

```bash
# From a Seurat object
Rscript tools/export_embedding.R \
  --input      /path/to/seurat.rds \
  --reductions pca,umap \
  --labels     guide,Phase \
  --scores     het_cv \
  --n          6000 \
  --source     "MDA-MB-231 CRISPRi Perturb-seq" \
  --caption    "6,000 cells from ..."

# ...or from an AnnData object
python3 tools/export_embedding.py --input adata.h5ad --labels guide,phase

# Always validate before committing
python3 tools/validate_embedding.py assets/data/embedding.json
```

Both exporters emit **only** reduction coordinates and the metadata columns you
name explicitly — never the counts matrix, barcodes, or unnamed columns. See
`assets/data/SCHEMA.md` for the format and for what must never go in the file.

> Use **cell-line** data. Individual-level TCGA expression is dbGaP-controlled
> and a public web page is not a permitted destination for it.

## Regenerating images

```bash
python3 tools/make_placeholder.py          # synthetic embedding.json
python3 tools/make_og_image.py             # assets/img/og.svg
# then rasterise og.svg -> og.png with headless Chrome (see that file's docstring)
```

## Deploying

The repo is a GitHub user site: pushing to `main` publishes the root of the repo
at `https://ssobt.github.io`. `.nojekyll` keeps Jekyll out of the way.

## Notes for future edits

- Header and footer markup is duplicated across the three HTML files on purpose.
  With no build step, duplication is more robust than a runtime include, and it
  keeps the pages working with JavaScript disabled.
- All colour, spacing, and type decisions live as custom properties at the top of
  `site.css`. Change them there, not inline.
- The categorical palette is Okabe–Ito, chosen to stay distinguishable under the
  common forms of colour vision deficiency. If you add series, extend
  `CATEGORICAL` in `scatter.js` rather than hard-coding colours.
