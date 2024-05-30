#!/usr/bin/env Rscript
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
# use for transpose of sparse matrix.
library(Matrix)
library(Seurat)
library(BPCells)
library(optparse)
library(tidyverse)

# * set up arguments
args <- OptionParser(usage = "Transform SnapATAC2 AnnData to Seurat v5 R object.",
  description = "

 NOTE:
 - R version >= 4.0
 - Seurat version >= 5.0
 - SnapATAC2 >= 2.4
 - cellmeta file is needed

 Usage examples
 0. Get help
    Rscript SnapATACAnnData2Seurat.R -h
 1. Directly transform without downsampling:
    Rscript SnapATACAnnData2Seurat.R -f [input.h5ad] -o [outfnm]
 2. Using BPCells in Seurat5
    - seurat object: [outfnm] and _mat subdirectory for BPCells
    Rscript SnapATACAnnData2Seurat.R -f [input.h5ad] -o [outfnm] --useBPCells
 3. Set python conda env to have SnapATAC2
    Rscript SnapATACAnnData2Seurat.R -f [input.h5ad] -o [outfnm] \
                             --conda [conda_path/bin/conda] --condaenv [env]
 4. Downsampling based on a column from a cell meta file
    Rscript SnapATACAnnData2Seurat.R -f [input.h5ad] -o [outfnm] --downsample \
                             -k [barcode_col] --dscol [cluster_col] --nds 50 \
                             -m [cellmetafile]
 5. Downsampling to a given number of cells.
    Rscript SnapATACAnnData2Seurat.R -f [input.h5ad] -o [outfnm] --downsample \
                             -k barcode -n 40000"
) |>
  add_option(opt_str = c("-f", "--sa2fnm"),
    type = "character", help = "snapatac2 AnnData file") |>
  add_option(opt_str = c("-o", "--outfnm"),
    type = "character", help = "seurat object output file name") |>
  add_option(opt_str = c("--conda"),
    type = "character", help = "conda binary path, default NULL") |>
  add_option(opt_str = c("--condaenv"),
    type = "character", default = NULL,
    help = "conda env for snapatac2, default NULL") |>
  add_option(opt_str = c("-m", "--sa2meta"),
    type = "character", default = NULL,
    help = "meta data file, csv format, default NULL") |>
  add_option(opt_str = c("--useBPCells"),
    type = "logical", default = FALSE,
    action = "store_true",
    help = "if using BPCells for large mat, default FALSE") |>
  add_option(opt_str = c("-d", "--downsample"),
    type = "logical", default = FALSE, action = "store_true",
    help = "if need to downsampling before transform, default FALSE") |>
  add_option(opt_str = c("-k", "--keycol"),
    type = "character", default = "barcode",
    help = "column of barcode in sa2meta, default barcode") |>
  add_option(opt_str = c("--dscol"),
    type = "character", default = NULL,
    help = "on which sa2meta column for downsample, default NULL") |>
  add_option(opt_str = c("--nds"), type = "integer",
    default = 50, help = "# of cells per group for downsample, defaualt 50") |>
  add_option(opt_str = c("-n", "--nmax"), type = "integer",
    default = 200000, help = "# cells at most, default 200000") |>
  parse_args(object = _)

## # * DEBUG
## datadir <- file.path("/tscc/projects/ps-renlab2/szu/projects",
##   "amb_pairedtag", "03.integration", "src/test/resource")
## args$sa2fnm <- file.path(datadir, "sa2ann_test.h5ad")
## args$outfnm <- file.path(datadir, "sa2ann_test.seu.rds")
## args$sa2meta <- file.path(datadir, "sa2ann_cellmeta.csv")
## args$keycol <- "barcode"
## args$dscol <- "brainregion"
## args$nds <- 50
## args$nmax <- 5000
## args$downsample <- FALSE
## args$conda <- "~/miniforge3/bin/mamba"
## args$condaenv <- "sa2"
## args$useBPCells <- FALSE

# * functions
getBarcodeFromAnn <- function(ann, backed = 'r') {
  if (backed == 'r') {
    barcode_ann <- ann$obs_names
  } else {
    barcode_ann <- ann$obs_names$to_list()
  }
  return(barcode_ann)
}

getFeatureFromAnn <- function(ann, backed = 'r') {
  if (backed == 'r') {
    r <- ann$var_names
  } else {
    r <- ann$var_names$to_list()
  }
  return(r)
}

checkMeta <- function(ann, meta, backed = 'r') {
  r <- meta
  barcode_ann <- getBarcodeFromAnn(ann, backed)
  if(nrow(meta) > length(barcode_ann)) {
    message("More barcodes in meta than in ann.")
    message("Only consider the ones in ann.")
    if (!all(barcode_ann %in% rownames(meta))) {
      stop("Some barcodes are not in the meta.")
    }
    r <- meta[barcode_ann, ]
  } else {
    if(!all(rownames(meta) %in% barcode_ann)) {
      stop("Some barcodes in meta are not in the ann.")
    }
  }
  return(r)
}

tos5 <- function(ann,
                 outfnm,
                 backed = 'r',
                meta = NULL,
                saveAsBPCells = FALSE) {
  barcode_ann <- getBarcodeFromAnn(ann, backed)
  features_ann <- getFeatureFromAnn(ann, backed)
  if (!is.null(meta)) {
    message("Cellmeta is not NULL, will save is to seurat too.")
    if (!all(rownames(meta) == barcode_ann)) {
      stop("meta and seurat barcodes are not matched.")
    }
  }
  outdir <- dirname(outfnm)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  message("Transform snapatac2 anndaata to Seurat",
    " with MatrixExtra package.")
  message("Treat X from ann as count, and force it to integer.")
  # ann and seurat rows are different.
  mat <- t(ann$X)
  rownames(mat) <- features_ann
  colnames(mat) <- barcode_ann
  mat <- MatrixExtra::as.csc.matrix(x = mat) |>
    BPCells::convert_matrix_type(matrix = _, type = "uint32_t")
  if (saveAsBPCells) {
    message("counts will be saved under ",
      file.path(outdir, "_mat"), " using BPCells.")
    outs5matdir <- file.path(outdir, "_mat")
    dir.create(outs5matdir, showWarnings = FALSE)
    BPCells::write_matrix_dir(mat = mat, outs5matdir,
      overwrite = TRUE)
    mat <- BPCells::open_matrix_dir(outs5matdir)
  } else {
    mat <- as(object = mat, Class = "dgCMatrix")
  }
  ann_meta <- ann$obs |>
    x => `rownames<-`(x, colnames(mat))
  ann_meta$barcode <- rownames(ann_meta)
  if (!is.null(meta)) {
    ann_meta <- merge(ann_meta, meta, by = "barcode") |>
      x => `rownames<-`(x, x$barcode) |>
      x => x[colnames(mat), ]
  }
  s5 <- Seurat::CreateSeuratObject(counts = mat, meta.data = ann_meta)
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
sa2 <- reticulate::import(module = "snapatac2")

# meta data from outside instead of from SnapATAC2
message("load SnapATAC2 meta csv file: ", args$sa2meta,
  " with data.table.")
meta <- data.table::fread(
  file = args$sa2meta, sep = ",",
  header = TRUE, data.table = FALSE) |>
  x => `rownames<-`(x, x[[args$keycol]])
if (args$keycol != "barcode") {
  message("keycol from meta is not barcode.")
  message("Then add barcode column in the meta.")
  meta$barcode <- meta[[args$keycol]]
}

if (!args$downsample) {
  message("Transform ann to seurat without dowmsample")
  message("load SnapATAC2 anndata file: ",
    args$sa2fnm, " into memory.")
  raw_ann <- sa2$read(filename = args$sa2fnm, backed = NULL)
  meta <- checkMeta(raw_ann, meta, backed = 'm')
  tos5(ann = raw_ann, outfnm = args$outfnm, meta = meta,
    saveAsBPCells = args$useBPCells, backed = 'm')
  message("Quit the script. Good Luck!")
  quit(save = "no", status = 0, )
}

message("load SnapATAC2 anndata file: ",
  args$sa2fnm, " in read-only mode.")
raw_ann <- sa2$read(filename = args$sa2fnm, backed = "r")
meta <- checkMeta(raw_ann, meta, backed = 'r')

if (is.null(args$dscol) || is.null(meta)) {
  message("Downsample before transform.")
  message(str_glue("Downsample to {args$nmax} barcodes."))
  barcodes <- sample(getBarcodeFromAnn(raw_ann, backed = 'r'),
    size = args$nmax, replace = FALSE)
} else {
  message(str_glue("Downsample based on {args$dscol} col in meta."))
  message(str_glue("Per group: {args$nds} at most."))
  barcodes <- dplyr::group_by(.data = meta, dplyr::across(args$dscol)) |>
    dplyr::slice_sample(n = args$nds) |>
    x => x[[args$keycol]]
}

message(str_glue("{length(barcodes)} barcodes are got."))
local({
  tmpf <- withr::local_tempfile(
    tmpdir = tempdir(),
    fileext = ".h5ad"
  )
  message(str_glue("save subset of ann to tmp file: {tmpf}."))
  message("The tmp file will be deleted at end.")
  sub_ann <- raw_ann$subset(obs_indices = barcodes,
    out = tmpf)
  sub_ann <- sub_ann$to_memory()
  sub_meta <- meta[barcodes, ]
  tos5(ann = sub_ann, outfnm = args$outfnm, meta = sub_meta,
    saveAsBPCells = args$useBPCells, backed = 'm')
})
message("Quit the script. Good Luck!")
