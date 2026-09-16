#!/usr/bin/env Rscript
# =============================================================================
# export_cohort.R — expression matrix  ->  assets/data/cohort.json
# -----------------------------------------------------------------------------
# Reproduces the patient-cohort analysis in Figure 1C-D of Woo & Sobti et al.
# (2026) and writes it in the format the site's cohort figure reads:
#
#   1. PCA over samples.
#   2. For each stratifier gene, split samples into bottom (q25) and top (q75)
#      expression quartiles.
#   3. For each sample, distance to ITS OWN GROUP'S centroid across the top
#      `--pcs` principal components.
#
# Emits ONLY PC coordinates, per-gene expression, quartile labels and
# distances. Sample identifiers are never written — patients are exported as
# anonymous indices.
#
# TCGA note: gene expression quantification from the GDC is OPEN ACCESS. The
# controlled tier is raw sequence and germline variants, neither of which is
# involved here. Publishing PCA coordinates derived from TCGA expression is
# what the paper itself does. Do not add barcodes to the output.
#
# Usage
#   Rscript tools/export_cohort.R \
#     --input      tcga_brca_fpkm.rds \
#     --genes      RNF8,MIS18A \
#     --controls   ADORA2B,ARHGEF5 \
#     --output     assets/data/cohort.json \
#     --pcs        50 \
#     --log        true \
#     --cohort     "TCGA-BRCA" \
#     --caption    "1,200 breast tumours from TCGA ..."
#
# --input is an .rds holding a numeric matrix with GENES IN ROWS and SAMPLES
# IN COLUMNS (rownames = gene symbols). Only --input and --genes are required.
# =============================================================================

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required:  install.packages('jsonlite')", call. = FALSE)
}

parse_args <- function(argv) {
  d <- list(input = NULL, output = "assets/data/cohort.json", genes = NULL,
            controls = "", pcs = "50", topvar = "2000", log = "true",
            cohort = "", caption = "", digits = "3")
  i <- 1
  while (i <= length(argv)) {
    k <- sub("^--", "", argv[[i]])
    if (!k %in% names(d)) stop("unknown argument: ", argv[[i]], call. = FALSE)
    if (i + 1 > length(argv)) stop("missing value for --", k, call. = FALSE)
    d[[k]] <- argv[[i + 1]]; i <- i + 2
  }
  if (is.null(d$input)) stop("--input <matrix.rds> is required", call. = FALSE)
  if (is.null(d$genes)) stop("--genes <SYMBOL,SYMBOL> is required", call. = FALSE)
  d$pcs <- as.integer(d$pcs); d$topvar <- as.integer(d$topvar)
  d$digits <- as.integer(d$digits); d$log <- tolower(d$log) %in% c("true","1","yes")
  d
}

split_csv <- function(x) if (!nzchar(x)) character(0) else trimws(strsplit(x, ",")[[1]])

args <- parse_args(commandArgs(trailingOnly = TRUE))

# --- load --------------------------------------------------------------------

message("Reading ", args$input, " ...")
mat <- readRDS(args$input)
if (!is.matrix(mat)) mat <- as.matrix(mat)
if (is.null(rownames(mat))) stop("matrix needs gene symbols as rownames", call. = FALSE)

candidates <- split_csv(args$genes)
controls   <- split_csv(args$controls)
all_genes  <- c(candidates, controls)

missing <- setdiff(all_genes, rownames(mat))
if (length(missing)) {
  stop("gene(s) not in the matrix: ", paste(missing, collapse = ", "), call. = FALSE)
}

n <- ncol(mat)
message("Cohort: ", n, " samples x ", nrow(mat), " genes.")
if (n < 100) warning("fewer than 100 samples; quartiles will be very small.")

# --- PCA ---------------------------------------------------------------------
# Log-transform, restrict to the most variable genes, then center and scale —
# the usual bulk-expression PCA recipe.

expr_mat <- if (args$log) log2(mat + 1) else mat

vars <- apply(expr_mat, 1, stats::var)
keep <- head(order(vars, decreasing = TRUE), min(args$topvar, nrow(expr_mat)))
sub  <- expr_mat[keep, , drop = FALSE]
sub  <- sub[apply(sub, 1, stats::var) > 0, , drop = FALSE]
message("PCA on the ", nrow(sub), " most variable genes.")

pca <- stats::prcomp(t(sub), center = TRUE, scale. = TRUE)
npc <- min(args$pcs, ncol(pca$x))
scores <- pca$x[, seq_len(npc), drop = FALSE]
varexp <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(npc)]
message(sprintf("Top %d PCs retained (PC1 %.1f%%, PC2 %.1f%% of variance).",
                npc, 100 * varexp[1], 100 * varexp[2]))

# --- per-gene stratification --------------------------------------------------

Q25 <- 0L; MID <- 1L; Q75 <- 2L

genes_out <- list()
for (g in all_genes) {
  e <- as.numeric(mat[g, ])
  cuts <- stats::quantile(e, c(0.25, 0.75), na.rm = TRUE)

  quart <- rep(MID, n)
  quart[e <= cuts[1]] <- Q25
  quart[e >= cuts[2]] <- Q75

  dist <- numeric(n)
  for (q in c(Q25, MID, Q75)) {
    idx <- which(quart == q)
    if (!length(idx)) next
    # Distance to this group's own centroid, across the retained PCs.
    cen <- colMeans(scores[idx, , drop = FALSE])
    dist[idx] <- sqrt(rowSums(sweep(scores[idx, , drop = FALSE], 2, cen, "-")^2))
  }

  med25 <- stats::median(dist[quart == Q25])
  med75 <- stats::median(dist[quart == Q75])
  w <- suppressWarnings(
    stats::wilcox.test(dist[quart == Q75], dist[quart == Q25])$p.value)

  message(sprintf("  %-10s median q25=%7.2f  q75=%7.2f  ratio=%.2f  P=%.3g",
                  g, med25, med75, med75 / med25, w))

  genes_out[[g]] <- list(
    label      = g,
    role       = if (g %in% controls) "control" else "candidate",
    expression = round(e, args$digits),
    quartile   = as.integer(quart),
    distance   = round(dist, 2),
    published  = list(centroidP = formatC(w, format = "g", digits = 3))
  )
}

# --- assemble ----------------------------------------------------------------
# Interleave PC1/PC2 as [x0,y0,x1,y1,...]; no sample names are emitted.

coords <- as.numeric(round(t(scores[, 1:2, drop = FALSE]), args$digits))

caption <- if (nzchar(args$caption)) args$caption else {
  paste0(n, " tumours", if (nzchar(args$cohort)) paste0(" from ", args$cohort) else "",
         ". Patients split into bottom (q25) and top (q75) expression quartiles of the ",
         "selected gene; spread measured as distance to the group centroid across the ",
         "top ", npc, " principal components.")
}

payload <- list(
  meta = list(
    n = n, cohort = args$cohort, synthetic = FALSE, dims = npc,
    caption = caption,
    method = paste0("PCA over the ", nrow(sub), " most variable genes; distance to ",
                    "quartile-group centroid across the top ", npc, " PCs."),
    varianceExplained = round(varexp[1:2], 4),
    exportedBy = "tools/export_cohort.R"
  ),
  coords    = list(pc = coords),
  geneOrder = all_genes,
  genes     = genes_out
)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, args$output, auto_unbox = TRUE, digits = args$digits)

kb <- round(file.info(args$output)$size / 1024)
message("Wrote ", args$output, " (", kb, " KB, ", n, " patients).")
message("NOTE: `published` now holds P-values computed from THIS data. If you want ",
        "the figure to cite the preprint's numbers instead, edit them in the JSON.")
