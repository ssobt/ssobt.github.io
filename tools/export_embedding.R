#!/usr/bin/env Rscript
# =============================================================================
# export_embedding.R — Seurat object  ->  assets/data/embedding.json
# -----------------------------------------------------------------------------
# Emits ONLY reduction coordinates and the metadata columns you name.
# The counts matrix, barcodes, and every unnamed metadata column stay behind.
#
# Usage
#   Rscript tools/export_embedding.R \
#     --input      seurat.rds \
#     --output     assets/data/embedding.json \
#     --reductions pca,umap \
#     --labels     guide,Phase \
#     --scores     het_cv \
#     --n          6000 \
#     --source     "MDA-MB-231 CRISPRi Perturb-seq" \
#     --caption    "6,000 cells from ..."
#
# Only --input is required. Defaults are listed in the parser below.
# =============================================================================

suppressPackageStartupMessages({
  ok <- requireNamespace("jsonlite", quietly = TRUE)
})
if (!ok) stop("jsonlite is required:  install.packages('jsonlite')", call. = FALSE)

# --- argument parsing --------------------------------------------------------

parse_args <- function(argv) {
  defaults <- list(
    input      = NULL,
    output     = "assets/data/embedding.json",
    reductions = "pca,umap",
    cloud      = "",
    labels     = "",
    scores     = "",
    n          = "6000",
    seed       = "42",
    source     = "",
    caption    = "",
    digits     = "3"
  )
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[[i]])
    if (!key %in% names(defaults)) stop("unknown argument: ", argv[[i]], call. = FALSE)
    if (i + 1 > length(argv)) stop("missing value for --", key, call. = FALSE)
    defaults[[key]] <- argv[[i + 1]]
    i <- i + 2
  }
  if (is.null(defaults$input)) {
    stop("--input <seurat.rds> is required", call. = FALSE)
  }
  defaults$n      <- as.integer(defaults$n)
  defaults$seed   <- as.integer(defaults$seed)
  defaults$digits <- as.integer(defaults$digits)
  defaults
}

split_csv <- function(x) {
  if (!nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

# --- load --------------------------------------------------------------------

message("Reading ", args$input, " ...")
obj <- readRDS(args$input)

if (!inherits(obj, "Seurat")) {
  stop("--input must be a saved Seurat object (got class: ",
       paste(class(obj), collapse = "/"), ")", call. = FALSE)
}

reductions <- split_csv(args$reductions)
# The hero's 3D intro cloud: first three dims of this reduction. Defaults to
# the first listed reduction (normally PCA, where PC3 is meaningful).
cloud_red <- if (nzchar(args$cloud)) args$cloud else reductions[1]
label_cols <- split_csv(args$labels)
score_cols <- split_csv(args$scores)

available <- names(obj@reductions)
missing_red <- setdiff(reductions, available)
if (length(missing_red)) {
  stop("reduction(s) not found: ", paste(missing_red, collapse = ", "),
       "\n  available: ", paste(available, collapse = ", "), call. = FALSE)
}

meta <- obj@meta.data
missing_meta <- setdiff(c(label_cols, score_cols), colnames(meta))
if (length(missing_meta)) {
  stop("metadata column(s) not found: ", paste(missing_meta, collapse = ", "),
       "\n  available: ", paste(colnames(meta), collapse = ", "), call. = FALSE)
}

# --- deterministic subsample -------------------------------------------------

n_total <- ncol(obj)
set.seed(args$seed)
keep <- if (n_total > args$n) sort(sample.int(n_total, args$n)) else seq_len(n_total)
n <- length(keep)
message("Exporting ", n, " of ", n_total, " cells (seed ", args$seed, ").")

# --- coordinates -------------------------------------------------------------

round_flat <- function(m, digits) {
  # Interleave to [x0,y0,x1,y1,...] and round.
  as.numeric(round(t(m[, 1:2, drop = FALSE]), digits))
}

# 3D cloud for the hero intro — skipped if the reduction has <3 dimensions.
coords3d <- NULL
if (cloud_red %in% available) {
  cloud_emb <- obj@reductions[[cloud_red]]@cell.embeddings
  if (ncol(cloud_emb) >= 3) {
    coords3d <- as.numeric(round(t(cloud_emb[keep, 1:3, drop = FALSE]), args$digits))
    message("3D hero cloud from '", cloud_red, "' dims 1-3.")
  } else {
    message("NOTE: '", cloud_red, "' has <3 dimensions; hero will start flat.")
  }
}

coords <- list()
for (r in reductions) {
  emb <- obj@reductions[[r]]@cell.embeddings
  if (ncol(emb) < 2) stop("reduction '", r, "' has fewer than 2 dimensions", call. = FALSE)
  coords[[r]] <- round_flat(emb[keep, , drop = FALSE], args$digits)
}

# --- label sets (categorical) ------------------------------------------------

label_sets <- list()
for (col in label_cols) {
  v <- meta[[col]][keep]
  f <- if (is.factor(v)) droplevels(v) else factor(v)
  lv <- levels(f)
  if (length(lv) > 24) {
    warning("column '", col, "' has ", length(lv),
            " levels; the legend will be unwieldy. Consider collapsing it.")
  }
  label_sets[[col]] <- list(
    label  = col,
    levels = as.character(lv),
    values = as.integer(f) - 1L            # zero-indexed for JS
  )
}

# --- scores (continuous) -----------------------------------------------------

scores <- list()
for (col in score_cols) {
  v <- as.numeric(meta[[col]][keep])
  if (anyNA(v)) {
    warning("column '", col, "' contains NA; replacing with the column median.")
    v[is.na(v)] <- stats::median(v, na.rm = TRUE)
  }
  scores[[col]] <- list(
    label  = col,
    range  = round(as.numeric(range(v)), args$digits),
    values = round(v, args$digits)
  )
}

# --- assemble ----------------------------------------------------------------

caption <- if (nzchar(args$caption)) args$caption else {
  paste0(n, " cells. ", toupper(reductions[length(reductions)]),
         " embedding; color by ",
         if (length(label_cols)) label_cols[1] else "cluster", ".")
}

payload <- list(
  meta = list(
    n               = n,
    source          = args$source,
    caption         = caption,
    heroCaption     = paste0(format(n, big.mark = ","), " cells · ",
                             toupper(reductions[length(reductions)])),
    embeddingLabels = setNames(as.list(toupper(reductions)), reductions),
    exportedBy      = "tools/export_embedding.R",
    seed            = args$seed
  ),
  coords    = coords,
  coords3d  = coords3d,
  labelSets = label_sets,
  scores    = scores
)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, args$output, auto_unbox = TRUE, digits = args$digits)

size_kb <- round(file.info(args$output)$size / 1024)
message("Wrote ", args$output, " (", size_kb, " KB, ", n, " cells).")
if (size_kb > 1024) {
  message("NOTE: over 1 MB. Lower --n so the page still paints quickly.")
}
message("Validate with:  python3 tools/validate_embedding.py ", args$output)
