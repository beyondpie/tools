#!/usr/bin/env Rscript
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
# use for transpose of sparse matrix.
library(Matrix)
library(Seurat)
library(BPCells)
library(optparse)
library(tidyverse)

# * set up arguments
args <- OptionParser(usage = "Transform Scanpy AnnData to Seurat v5 R object.",
  description = "

 NOTE:
 - R version >= 4.0
 - Seurat version >= 5.0


 Usage examples
 0. Get help
    Rscript SnapATACAnnData2Seurat.R -h
 1. Directly transform without downsampling:
    Rscript ScanpyAnnData2Seurat.R -f [input.h5ad] -o [outfnm]
 2. Using BPCells in Seurat5
    - seurat object: [outfnm] and _mat subdirectory for BPCells
    Rscript ScanpyAnnData2Seurat.R -f [input.h5ad] -o [outfnm] --useBPCells
 3. Set python conda env to have anndata
    Rscript ScanpyAnnData2Seurat.R -f [input.h5ad] -o [outfnm] \
                             --conda [conda_path/bin/conda] --condaenv [env]
 4. Downsampling to a given number of cells.
    Rscript ScanpyAnnData2Seurat.R -f [input.h5ad] -o [outfnm] --downsample \
                             -k barcode -n 40000"
) |>
  add_option(opt_str = c("-f", "--annfnm"),
    type = "character", help = "snapatac2 AnnData file") |>
  add_option(opt_str = c("-o", "--outfnm"),
    type = "character", help = "seurat object output file name") |>
  add_option(opt_str = c("--conda"),
    type = "character", help = "conda binary path, default NULL") |>
  add_option(opt_str = c("--condaenv"),
    type = "character", default = NULL,
    help = "conda env for snapatac2, default NULL") |>
  add_option(opt_str = c("--useBPCells"),
    type = "logical", default = FALSE,
    action = "store_true",
    help = "if using BPCells for large mat, default FALSE") |>
  add_option(opt_str = c("-d", "--downsample"),
    type = "logical", default = FALSE, action = "store_true",
    help = "if need to downsampling before transform, default FALSE") |>
  add_option(opt_str = c("--dscol"),
    type = "character", default = NULL,
    help = "on which annmeta column for downsample, default NULL") |>
  add_option(opt_str = c("--nds"), type = "integer",
    default = 50, help = "# of cells per group for downsample, defaualt 50") |>
  add_option(opt_str = c("-n", "--nmax"), type = "integer",
    default = 200000, help = "# cells at most, default 200000") |>
  add_option(opt_str = c("--matGroup"), type = "character",
    default = "X",
    help = "mat field, usually X, but layers/rawcount in Allen") |>
  parse_args(object = _)

## # * DEBUG
## datadir <- file.path("/tscc/projects/ps-renlab2/szu/projects/",
##   "amb_pairedtag", "03.integration", "src/test/resource")
## args$annfnm <- file.path(datadir, "scanpy_ann.h5ad")
## args$outfnm <- file.path(datadir, "sc.ann2seurat.rds")
## args$dscol <- "sex"
## args$nds <- 50
## args$nmax <- 100
## args$downsample <- FALSE
## args$conda <- "/tscc/nfs/home/szu/miniforge3/bin/conda"
## args$condaenv <- "sa2"
## args$useBPCells <- FALSE
## args$matGroup <- "layers/rawcount"

# * functions
tos5 <- function(ann,
                 outfnm,
                 matGroup = "X",
                saveAsBPCells = FALSE) {
  barcode_ann <- ann$obs_names$to_list()
  features_ann <- ann$var_names$to_list()
  outdir <- dirname(outfnm)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  message("Transform snapatac2 anndaata to Seurat",
    " with MatrixExtra package.")
  message("Treat ", matGroup,
    " from ann as count, and force it to integer.")
  # ann and seurat rows are different.
  if(matGroup == "X") {
    mat <- t(ann$X)
  } else {
    mat <- t(ann$layers[[gsub("layers/", "", matGroup)]])
  }
  rownames(mat) <- features_ann
  colnames(mat) <- barcode_ann
  mat <- MatrixExtra::as.csc.matrix(x = mat) |>
    BPCells::convert_matrix_type(matrix = _, type = "uint32_t")
  if (saveAsBPCells) {
    message(
      "counts will be saved under ",
      file.path(outdir, "_mat"), " using BPCells."
    )
    outs5matdir <- file.path(outdir, "_mat")
    dir.create(outs5matdir, showWarnings = FALSE)
    BPCells::write_matrix_dir(
      mat = mat, outs5matdir,
      overwrite = TRUE
    )
    mat <- BPCells::open_matrix_dir(outs5matdir)
  } else {
    mat <- as(object = mat, Class = "dgCMatrix")
  }
  s5 <- Seurat::CreateSeuratObject(counts = mat,
    meta.data = ann$obs)
  message("Save Seurat to ", outfnm)
  saveRDS(s5, outfnm)
  message("Seurat object saved.")
}

# * main
if (!is.null(args$conda)) {
  message("Using reticulate to setup conda: ",
    args$conda)
  if (is.null(args$condaenv)) {
    message("No conda env is set, use default one.")
  }
  reticulate::use_condaenv(
    condaenv = args$condaenv,
    conda = args$conda,
    required = TRUE)
}
ad <- reticulate::import(module = "anndata")
message("load anndata file: ", args$annfnm)
raw_ann <- ad$read_h5ad(filename = args$annfnm)
raw_barcodes <- raw_ann$obs_names$to_list()
meta <- raw_ann$obs |>
  x => `rownames<-`(x, raw_barcodes)
meta$barcode <- rownames(meta)

if (!args$downsample) {
  message("Transform ann to seurat without dowmsample.")
  tos5(ann = raw_ann, outfnm = args$outfnm,
    matGroup = args$matGroup,
    saveAsBPCells = args$useBPCells)
  message("Quit the script. Good Luck!")
  quit(save = "no", status = 0, )
}

if (is.null(args$dscol)) {
  message("Downsample before transform.")
  message(str_glue("Downsample to {args$nmax} barcodes."))
  barcodes <- sample(raw_barcodes, size = args$nmax, replace = FALSE)
} else {
  message(str_glue("Downsample based on {args$dscol} col in meta."))
  message(str_glue("Per group: {args$nds} at most."))
  barcodes <- dplyr::group_by(.data = meta, dplyr::across(args$dscol)) |>
    dplyr::slice_sample(n = args$nds) |>
    x => x$barcode
}

message(str_glue("{length(barcodes)} barcodes are got."))
sub_ann <- raw_ann[raw_barcodes %in% barcodes]$copy()
tos5(ann = sub_ann, outfnm = args$outfnm,
  matGroup = args$matGroup,
  saveAsBPCells = args$useBPCells)
message("Quit the script. Good Luck!")
