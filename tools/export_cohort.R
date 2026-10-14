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
#   4. An EMPIRICAL NULL from `--null-genes` randomly chosen expressed genes,
#      put through the identical procedure.
#
# Step 4 is not optional. With ~275 samples per quartile, roughly 80% of
# arbitrary genes reach P < 0.05 on this metric — a raw p-value says almost
# nothing. What matters is where a gene falls against the null, which is the
# same calibration the paper uses to nominate candidates (empirical z-score).
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
# --input is either
#   * an .rds holding a numeric matrix, genes in rows, samples in columns
#     (rownames = gene symbols), or
#   * a cBioPortal-style tab-delimited matrix such as TCGA's
#     data_RNA_Seq_v2_expression_median.txt, whose first two columns are
#     Hugo_Symbol and Entrez_Gene_Id and whose remaining columns are samples.
# Only --input and --genes are required.
# =============================================================================

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required:  install.packages('jsonlite')", call. = FALSE)
}

parse_args <- function(argv) {
  d <- list(input = NULL, output = "assets/data/cohort.json", genes = NULL,
            controls = "", pcs = "50", topvar = "0", log = "true",
            cohort = "", caption = "", digits = "3", null_genes = "200",
    min_expr = "1")
  i <- 1
  while (i <= length(argv)) {
    k <- gsub("-", "_", sub("^--", "", argv[[i]]))   # --null-genes -> null_genes
    if (!k %in% names(d)) stop("unknown argument: ", argv[[i]], call. = FALSE)
    if (i + 1 > length(argv)) stop("missing value for --", k, call. = FALSE)
    d[[k]] <- argv[[i + 1]]; i <- i + 2
  }
  if (is.null(d$input)) stop("--input <matrix.rds> is required", call. = FALSE)
  if (is.null(d$genes)) stop("--genes <SYMBOL,SYMBOL> is required", call. = FALSE)
  d$pcs <- as.integer(d$pcs); d$topvar <- as.integer(d$topvar)
  d$null_genes <- as.integer(d$null_genes); d$min_expr <- as.numeric(d$min_expr)
  d$digits <- as.integer(d$digits); d$log <- tolower(d$log) %in% c("true","1","yes")
  d
}

split_csv <- function(x) if (!nzchar(x)) character(0) else trimws(strsplit(x, ",")[[1]])

args <- parse_args(commandArgs(trailingOnly = TRUE))

# --- load --------------------------------------------------------------------

message("Reading ", args$input, " ...")

read_matrix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    m <- readRDS(path)
    if (!is.matrix(m)) m <- as.matrix(m)
    return(m)
  }
  # Tab-delimited (cBioPortal / TCGA) export.
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("reading a .txt matrix needs data.table:  install.packages('data.table')",
         call. = FALSE)
  }
  dt <- data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  sym_col <- intersect(c("Hugo_Symbol", "gene", "Gene", "GENE_SYMBOL"), names(dt))[1]
  if (is.na(sym_col)) sym_col <- names(dt)[1]
  drop <- intersect(c(sym_col, "Entrez_Gene_Id", "EntrezID"), names(dt))
  syms <- as.character(dt[[sym_col]])

  m <- as.matrix(dt[, !drop, with = FALSE])
  storage.mode(m) <- "numeric"

  # Gene symbols repeat in these files (multiple Entrez IDs per symbol). Keep
  # the copy with the highest mean expression rather than an arbitrary one.
  ok <- !is.na(syms) & nzchar(syms)
  m <- m[ok, , drop = FALSE]; syms <- syms[ok]
  if (anyDuplicated(syms)) {
    ord <- order(syms, -rowMeans(m, na.rm = TRUE))
    m <- m[ord, , drop = FALSE]; syms <- syms[ord]
    keep <- !duplicated(syms)
    message("  collapsed ", sum(!keep), " duplicate gene symbols (kept highest-expressed)")
    m <- m[keep, , drop = FALSE]; syms <- syms[keep]
  }
  rownames(m) <- syms
  m
}

mat <- read_matrix(args$input)
if (is.null(rownames(mat))) stop("matrix needs gene symbols as rownames", call. = FALSE)

# Samples with all-NA rows break prcomp; drop them up front.
bad <- colSums(is.na(mat)) == nrow(mat)
if (any(bad)) { message("  dropping ", sum(bad), " all-NA samples"); mat <- mat[, !bad, drop = FALSE] }
mat[is.na(mat)] <- 0

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

# topvar = 0 means "use every gene". Restricting to highly variable genes is a
# defensible choice but it changes which axes dominate, so the default here is
# the whole transcriptome.
if (args$topvar > 0) {
  vars <- apply(expr_mat, 1, stats::var)
  keep <- head(order(vars, decreasing = TRUE), min(args$topvar, nrow(expr_mat)))
  sub  <- expr_mat[keep, , drop = FALSE]
} else {
  sub <- expr_mat
}
sub <- sub[apply(sub, 1, stats::var) > 0, , drop = FALSE]
message("PCA on ", nrow(sub), " genes",
        if (args$topvar > 0) " (most variable)" else " (all)", ".")

pca <- stats::prcomp(t(sub), center = TRUE, scale. = TRUE)
npc <- min(args$pcs, ncol(pca$x))
scores <- pca$x[, seq_len(npc), drop = FALSE]
varexp <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(npc)]
message(sprintf("Top %d PCs retained (PC1 %.1f%%, PC2 %.1f%% of variance).",
                npc, 100 * varexp[1], 100 * varexp[2]))

# --- per-gene stratification --------------------------------------------------

Q25 <- 0L; MID <- 1L; Q75 <- 2L

# Distances to each quartile group's own centroid, for one stratifier vector.
spread_for <- function(e) {
  cuts <- stats::quantile(e, c(0.25, 0.75), na.rm = TRUE)
  if (!is.finite(cuts[1]) || cuts[1] == cuts[2]) return(NULL)
  quart <- rep(MID, n)
  quart[e <= cuts[1]] <- Q25
  quart[e >= cuts[2]] <- Q75
  if (sum(quart == Q25) < 30 || sum(quart == Q75) < 30) return(NULL)
  dist <- numeric(n)
  for (q in c(Q25, MID, Q75)) {
    idx <- which(quart == q)
    if (!length(idx)) next
    cen <- colMeans(scores[idx, , drop = FALSE])
    dist[idx] <- sqrt(rowSums(sweep(scores[idx, , drop = FALSE], 2, cen, "-")^2))
  }
  m25 <- stats::median(dist[quart == Q25]); m75 <- stats::median(dist[quart == Q75])
  list(quart = quart, dist = dist, ratio = m75 / m25, med25 = m25, med75 = m75,
       p = suppressWarnings(
         stats::wilcox.test(dist[quart == Q75], dist[quart == Q25])$p.value))
}

# --- empirical null -----------------------------------------------------------
set.seed(42)
expressed <- which(rowMeans(expr_mat) > args$min_expr)
pool <- setdiff(expressed, match(all_genes, rownames(mat)))
draws <- sample(pool, min(args$null_genes, length(pool)))
message("Building empirical null from ", length(draws), " random expressed genes ...")
null_ratio <- c(); null_sig <- 0
for (i in draws) {
  r <- spread_for(as.numeric(mat[i, ]))
  if (is.null(r)) next
  null_ratio <- c(null_ratio, r$ratio)
  if (!is.na(r$p) && r$p < 0.05) null_sig <- null_sig + 1
}
null_mean <- mean(null_ratio); null_sd <- stats::sd(null_ratio)
message(sprintf("  null spread ratio: mean %.3f  sd %.3f  (%.0f%% of random genes reach P<0.05)",
                null_mean, null_sd, 100 * null_sig / length(null_ratio)))

# --- per-gene stratification --------------------------------------------------

genes_out <- list()
for (g in all_genes) {
  e <- as.numeric(mat[g, ])
  r <- spread_for(e)
  if (is.null(r)) { message("  ", g, ": not enough spread to stratify; skipped"); next }
  z <- (r$ratio - null_mean) / null_sd
  pct <- 100 * mean(null_ratio < r$ratio)

  message(sprintf("  %-10s median q25=%7.2f  q75=%7.2f  ratio=%.3f  z=%+.2f  pctile=%4.1f%%  P=%.3g",
                  g, r$med25, r$med75, r$ratio, z, pct, r$p))

  genes_out[[g]] <- list(
    label      = g,
    role       = if (g %in% controls) "control" else "candidate",
    expression = round(e, args$digits),
    quartile   = as.integer(r$quart),
    distance   = round(r$dist, 2),
    stats      = list(
      ratio      = round(r$ratio, 4),
      z          = round(z, 2),
      percentile = round(pct, 1),
      wilcoxonP  = formatC(r$p, format = "g", digits = 3)
    )
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
  nullModel = list(
    nGenes     = length(null_ratio),
    meanRatio  = round(null_mean, 4),
    sdRatio    = round(null_sd, 4),
    pctRandomSignificant = round(100 * null_sig / length(null_ratio), 1),
    note = paste0("Spread ratios from ", length(null_ratio), " randomly chosen expressed ",
                  "genes put through the identical procedure. At this cohort size most ",
                  "arbitrary genes reach nominal significance, so a gene's z-score ",
                  "against this null is the meaningful quantity, not its p-value.")
  ),
  coords    = list(pc = coords),
  geneOrder = all_genes,
  genes     = genes_out
)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, args$output, auto_unbox = TRUE, digits = args$digits)

kb <- round(file.info(args$output)$size / 1024)
message("Wrote ", args$output, " (", kb, " KB, ", n, " patients).")
message("NOTE: per-gene `stats` are computed from THIS data. Any published values ",
        "the figure cites alongside them must be added separately.")
