#!/usr/bin/env Rscript

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
options(future.globals.maxSize = 10e9)
library(optparse)
library(Matrix)
library(Seurat)
library(tidyverse)
library(BPCells)
library(future)

# * set up arguments
args <- OptionParser(
  usage = "Seurat Transfer Label",
  description = "
 NOTE:
 - R version >= 4.2
 - Seurat version >= 5.0
 - Methods: as listed in Seurat::FindTransferAnchors

 Usage examples
 0. Get help
    Rscript TransferLabel.R -h
"
) |>
  add_option(
    opt_str = c("-r", "--ref"), type = "character",
    help = "reference seurat object fnm"
  ) |>
  add_option(
    opt_str = c("-q", "--query"), type = "character",
    help = "query seurat object fnm"
  ) |>
  add_option(
    opt_str = c("-o", "--outdir"), type = "character",
    help = "the directory to save the results"
  ) |>
  add_option(
    opt_str = c("-f", "--feature"), type = "character",
    default = NULL, help = "gene list file as feature"
  ) |>
  add_option(
    opt_str = c("--nvar"), type = "integer",
    default = 2000, help = " # of variable features"
  ) |>
  add_option(
    opt_str = c("--useref"), type = "logical",
    action = "store_ture", default = FALSE,
    help = "use variable features from ref"
  ) |>
  add_option(
    opt_str = c("--usequery"), type = "logical",
    action = "store_true", default = FALSE
  ) |>
  add_option(
    opt_str = c("--saveanchor"), type = "logical",
    action = "store_true", default = FALSE
  ) |>
  add_option(
    opt_str = c("--rerun"), type = "logical",
    action = "store_true", default = FALSE
  ) |>
  add_option(
    opt_str = c("-m", "--method"), type = "character",
    default = "cca",
    help = "method for transfer label"
  ) |>
  add_option(
    opt_str = c("-k", "--kanchor"), type = "integer",
    default = 5,
    help = "kanchor for finding anchors"
  ) |>
  add_option(
    opt_str = c("--threads"), type = "integer",
    default = 4,
    help = "threads to use"
  ) |>
  add_option(
    opt_str = c("-p", "--npca"), type = "integer",
    default = 50, help = "# of components for finding anchors"
  ) |>
  add_option(
    opt_str = c("-t", "--tfcol"), type = "character",
    help = "on which col of ref we transfer labels to query."
  ) |>
  parse_args(object = _)

# * DEBUG
# projdir <- "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
## projdir <- "/projects/ps-renlab2/szu/projects/amb_pairedtag"
## datadir <- file.path(
##   projdir, "03.integration", "src/test/resource"
## )
## args$ref <- file.path(datadir, "allenseu_ds.rds")
## args$query <- file.path(datadir, "ptseu_ds.rds")
## args$tfcol <- "cl"
## args$outdir <- datadir
## args$useref <- TRUE
## args$saveanchor <- TRUE
## args$kanchor <- 50
## args$method <- "pcaproject"
## args$feature <- file.path(
##   projdir, "meta", "AIT21_k8_markers.txt"
## )

message("threads: ", args$threads)
plan("multicore", workers = args$threads)

# * set up inputs
getNoZeroCountGene <- function(s5) {
  m <- s5@assays$RNA@layers$counts
  rownames(s5)[rowSums(m) >= 1]
}

message("Load ref seurat.")
ref <- readRDS(args$ref)
genes_ref_var <- getNoZeroCountGene(ref)

if (is.null(args$tfcol) || (!args$tfcol %in% colnames(ref@meta.data))) {
  stop("error: tfcol is NULL or not exist in ref.")
}

message("Load query seruat.")
query <- readRDS(args$query)
genes_query_var <- getNoZeroCountGene(query)

genes_all <- intersect(genes_ref_var, genes_query_var)
message(str_glue("{length(genes_all)} genes from both."))
ref <- ref[genes_all, ] |>
  NormalizeData(object = _)
query <- query[genes_all, ] |>
  NormalizeData(object = _)

message("prepare outdir.")
dir.create(args$outdir, showWarnings = FALSE)

message("Set up features for transfer label.")
if (!is.null(args$feature)) {
  message(str_glue("load features from file: {args$feature}."))
  features <- read.table(
    file = args$feature, header = FALSE
  )$V1
  features <- intersect(features, genes_all)
} else {
  message("Use variable features.")
  nvar <- args$nvar
  message(str_glue("Finding {nvar} variable features."))
  ref <- FindVariableFeatures(
    object = ref,
    selection.method = "vst",
    nfeatures = nvar
  )
  fea_ref <- VariableFeatures(ref)
  query <- FindVariableFeatures(
    object = query,
    selection.method = "vst",
    nfeatures = nvar
  )
  fea_query <- VariableFeatures(query)
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
    features <- union(fea_ref, fea_query)
  }
}
VariableFeatures(ref) <- features
VariableFeatures(query) <- features

# * set up anchor process
if (args$saveanchor) {
  anchorfnm <- file.path(args$outdir,
    str_glue("tf.anchors.with-{args$method}-kac{args$kanchor}.rds"))
}
if ((!args$rerun) && (file.exists(anchorfnm))) {
  message("Anchor file exist, and will load it.")
  anchors <- readRDS(anchorfnm)
} else {
  message(str_glue("FindTransferAnchors with {args$method} method."))
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
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
query <- TransferData(
  anchorset = anchors,
  refdata = ref@meta.data[[args$tfcol]],
  reference = NULL,
  query = query,
  weight.reduction = args$method,
  l2.norm = FALSE,
  dims = seq_len(args$npca),
  # default
  k.weight = 50,
  verbose = TRUE,
)
message("Transferred labels are written into query.")

queryfnm <- file.path(
  args$outdir,
  str_glue("query.with.tf-{args$method}-kac{args$kanchor}_on-{args$tfcol}.rds")
)
message(
  "Save query object with tflabels to: ",
  queryfnm
)
saveRDS(object = query, queryfnm)
