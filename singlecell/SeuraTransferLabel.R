#!/usr/bin/env Rscript

library(Matrix)
options(future.globals.maxSize = 5e9)
library(ggplot2)

import::from(.from = optparse, make_option,
  OptionParser, add_option, parse_args)
import::from(.from = stringr, str_glue)
import::from(.from = dplyr, group_by, slice_sample)
import::from(.from = Seurat, VariableFeatures,
  `VariableFeatures<-`,
  FindVariableFeatures,
  FindTransferAnchors,
  TransferData, DefaultAssay)

# This script is to provide a unified transfer label module
# - Only works for two dataset at one time
# - cca, rpca
# - support co-embedding plot
# - automatically label group based on votes
# - support consensus matrix plot
# - input format is Seurat

# NOTE:
# - R version >= 4.2
# - Seurat version >= 5.0

# * set up arguments
args <- OptionParser() |>
  add_option(opt_str = c("-r", "--ref"), type = "character",
    help = "reference seurat object fnm") |>
  add_option(opt_str = c("-q", "--query"), type = "character",
    help = "query seurat object fnm") |>
  add_option(opt_str = c("-o", "--outdir"), type = "character",
    help = "the directory to save the results") |>
  add_option(opt_str = c("--downsample"), type = "logical",
    action = "store_true", default = FALSE) |>
  add_option(opt_str = c("--dcolref"), type = "character",
    help = "col of ref for downsample") |>
  add_option(opt_str = c("--dcolquery"), type = "character",
    help = "col of query for downsample") |>
  add_option(opt_str = c("--ndref"), type = "integer",
    default = NULL,
    help = "# of cells per cluster for ref if downsampling") |>
  add_option(opt_str = c("--ndquery"), type = "integer",
    default = NULL,
    help = "# of cells per cluster for query if downsampling") |>
  add_option(opt_str = c("--nref"), type = "integer",
    default = 200000,
    help = "# of cells in total for ref") |>
  add_option(opt_str = c("--nquery"), type = "integer",
    default = 200000,
    help = "# of cells in total for query") |>
  add_option(opt_str = c("-f", "--feature"), type = "character",
    default = NULL, help = "gene list file as feature") |>
  add_option(opt_str = c("--nvar"), type = "integer",
    default = 2000, help = " # of variable features") |>
  add_option(opt_str = c("--useref"), type = "logical",
    action = "store_ture", default = FALSE,
    help = "use variable features from ref") |>
  add_option(opt_str = c("--usequery"), type = "logical",
    action = "store_true", default = FALSE) |>
  add_option(opt_str = c("--saveanchor"), type = "logical",
    action = "store_true", default = FALSE) |>
  add_option(opt_str = c("-m", "--method"), type = "character",
    default = "cca",
    help = "method for transfer label") |>
  add_option(opt_str = c("-k", "--kanchor"), type = "integer",
    defautl = 5,
    help = "kanchor for finding anchors") |>
  add_option(opt_str = c("-p", "--npca"), type = "integer",
    default = 50, help = "# of components for finding anchors") |>
  add_option(opt_str = c("-t", "--tfcol"), type = "character",
    help = "on which col of ref we transfer labels to query.") |>
    add_option(opt_str = c("--ordcolref"), type = "character",
      help = "ordered cols of ref for consensus matrix plot.") |>
    add_option(opt_str = c("--ordcolquery"), type = "character",
      help = "ordered cols of query for consensus matrix plot.") |>
    add_option(opt_str = c("--ecref"), type = "character",
      default = NULL,
      help = "extra col in ref for consensus matrix plot.") |>
    add_option(opt_str = c("--ecquery"), type = "character",
      default = NULL,
      help = "extra col in query for consensus matrix plot.") |>
    add_option(opt_str = c("--ordecref"), type = "character",
      default = NULL,
      help = "ordered extra cols of ref.") |>
    add_option(opt_str = c("--ordecquery"), type = "character",
      default = NULL,
      help = "ordered extra cols of query.") |>
    add_option(opt_str = c("--pltcoembed"), type = "logical",
      action = "store_true", default = FALSE,
      help = "if plot coembedding on two datasets.") |>
  parse_args(object = _)

# * DEBUG
## datadir <- file.path("~/git-recipes/mouseBrainAtlas",
##   "amb_pairedtag", "03.integration", "src/test/resource")
## args$ref <- file.path(datadir, "allen_ref.rds")
## args$query <- file.path(datadir, "pt_query.rds")
## args$tfcol <- "cl"
## args$outdir <- datadir
## args$downsample <- TRUE
## args$dcolref <- "cl"
## args$dcolquery <- "L3"
## args$useref <- TRUE
## args$saveanchor <- TRUE
## args$kanchor <- 50

# * set up inputs
eps <- 1e-6
# pretty slow
getVarGenes <- function(s5) {
  m <- s5@assays$RNA@layers$counts
  vapply(seq_len(nrow(m)), \(i) {
    (max(m[i, ]) - min(m[i, ])) > eps
  }, TRUE) |>
    x => rownames(s5)[x]
}

getNoZeroCountGene <- function(s5) {
  m <- s5@assays$RNA@layers$counts
  rownames(s5)[rowSums(m) > eps]
}

message("Run Seurat Transfer Label.")
ref <- readRDS(args$ref)
genes_ref_var <- getNoZeroCountGene(ref)

if (is.null(args$tfcol)) {
  stop("error: transfer col is NULL.")
}
if (!args$tfcol %in% colnames(ref@meta.data)) {
  stop(str_glue("error: {args$tfcol} is not in ref meta data."))
}

if ((!is.null(args$extrpltcol)) && (
  !args$extrpltcol %in% colnames(ref@meta.data))) {
  stop(str_glue("error: {args$extrpltcol} is not in ref meta data."))
}

query <- readRDS(args$query)
genes_query_var <- getNoZeroCountGene(query)

genes_all <- intersect(genes_ref_var, genes_query_var)
message(str_glue("{length(genes_all)} genes from both."))
ref <- ref[genes_all, ] |>
  Seurat::NormalizeData(object = _)
query <- query[genes_all, ] |>
  Seurat::NormalizeData(object = _)
rqlist <- list(ref = ref, query = query)
rm(ref, query)

# * outdir prepare
dir.create(args$outdir, showWarnings = FALSE)
message("Results will be put under: ", args$outdir)

# * downsampling

if (args$downsample) {
  rqlist2 <- lapply(names(rqlist), \(nm) {
    message(str_glue("Perform downsampling on {nm}."))
    s <- rqlist[[nm]]
    clusters <- s@meta.data[[
      args[[paste0("dcol", nm)]]
    ]]
    ntotal <- ncol(s)
    ncl <- length(table(clusters))
    nallow <- min(ntotal, args[[paste0("n", nm)]])
    nd <- args[[paste0("nd", nm)]] |>
      x => ifelse(test = is.null(x),
        floor(nallow / ncl), x)
    message(str_glue("In total, {ntotal} cells and {ncl} d-clusters."))
    message(str_glue("Downsample {nd} per cluster",
      " and {nallow} cells at most."))
    barcodes <- data.frame(
      barcode = colnames(s),
      cluster = clusters) |>
      group_by(cluster) |>
      slice_sample(n = nd) |>
      x => x$barcode
    message(str_glue("{length(barcodes)} barcodes are kept."))
    s[, barcodes]
  })
} else {
  rqlist2 <- rqlist
}
names(rqlist2) <- names(rqlist)
rm(rqlist)

# * set up features
if (!is.null(args$feature)) {
  message(str_glue("load features from file: {option$feature}."))
  features <- read.table(file = args$feature, header = FALSE,
    row.names = FALSE, col.names = FALSE)$V1
  features <- intersect(features, genes_all)
  VariableFeatures(rqlist2$ref) <- features
  VariableFeatures(rqlist2$query) <- features
} else {
  nvar <- args$nvar
  message(str_glue("Finding {nvar} variable features."))
  rqlist2$ref <- FindVariableFeatures(
    object = rqlist2$ref,
    selection.method = "vst",
    nfeatures = nvar
  )
  fea_ref <- VariableFeatures(rqlist2$ref)
  rqlist2$query <- FindVariableFeatures(
    object = rqlist2$query,
    selection.method = "vst",
    nfeatures = nvar
  )
  fea_query <- VariableFeatures(rqlist2$query)
  if (args$useref && (!args$usequery)) {
    message("use features from REF Variable features.")
    features <- fea_ref
  }
  if (args$useref && args$usequery) {
    message("use features from BOTH Variable features.")
    features <- intersect(fea_ref, fea_query)
    message(str_glue("{length(features)} features used."))
  }
  if ((!args$useref) && args$usequery) {
    message("use features from QUERY Variable features.")
    features <- use_query
  }
  if ((!args$useref) && (!args$usequery)) {
    message("use features from BOTH features.")
    features <- genes_all
  }
}
VariableFeatures(rqlist2$ref) <- features
VariableFeatures(rqlist2$query) <- features

# * set up anchor process
if (args$saveanchor) {
  anchorfnm <- file.path(args$outdir, "tf.anchors.rds")
}
if (file.exists(anchorfnm)) {
  message("Anchor file exist, and will load it.")
  anchors <- readRDS(anchorfnm)
} else {
  message(str_glue("FindTransferAnchors with {args$method} method."))
  anchors <- FindTransferAnchors(
    reference = rqlist2$ref,
    query = rqlist2$query,
    scale = TRUE,
    reduction = args$method,
    normalization.method = "LogNormalize",
    npcs = args$npca,
    dims = seq_len(args$npca),
    l2.norm = TRUE,
    k.anchor = args$kanchor,
    features = features,
    verbose = TRUE
  )
  if (args$saveanchor) {
    saveRDS(object = anchors, file = anchorfnm)
  }
}

# * set up transfer labels and annotations
message("Then transfer labels in single-cell level.")
rqlist2$query <- TransferData(
  anchorset = anchors,
  refdata = rqlist2$ref@meta.data[[args$tfcol]],
  reference = NULL,
  query = rqlist2$query,
  weight.reduction = ifelse(
    args$method != "cca", "pcaproject", "cca"),
  l2.norm = FALSE,
  dims = seq_len(args$npca),
  # default
  k.weight = 50,
  verbose = TRUE,
  )
message("Transferred labels are written into query.")

queryfnm <- file.path(args$outdir,
  str_glue("query.with.tflabels.{args$tfcol}.rds"))
message("Save query object with tflabels to: ",
  queryfnm)
saveRDS(object = rqlist2$query, queryfnm)

# * generate consensus matrix

# * plot consensus matrix
# * plot co-embedding
