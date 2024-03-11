#!/usr/bin/env Rscript
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

import::from(.from = optparse, add_option,
  OptionParser, add_option, parse_args)
import::from(.from = stringr, str_glue)

# Transform SnapATAC2 AnnObject to Seurat v5 object.
# NOTE:
# - R version >= 4.2
# - Seurat version >= 5.0
# - python version >= 3.8
# - Python package: snapatac2

# * set up arguments
args <- OptionParser() |>
  add_option(opt_str = c("-f", "--sa2fnm"),
    type = "character", help = "snapatac2 AnnData file") |>
  add_option(opt_str = c("-o", "--outdir"),
    type = "character", help = "seurat outdir") |>
  add_option(opt_str = c("--conda"),
    type = "character", help = "conda binary path") |>
  add_option(opt_str = c("--condaenv"),
    type = "character", default = NULL,
    help = "conda env for snapatac2") |>
  add_option(opt_str = c("-m", "--sa2meta"),
    type = "character", default = NULL,
    help = "meta data file, csv format") |>
  add_option(opt_str = c("--useBPCells"),
    type = "logical", default = FALSE,
    action = "store_true",
    help = "if using BPCells for large mat") |>
  add_option(opt_str = c("-d", "--downsample"),
    type = "logical", default = FALSE, action = "store_true",
    help = "if need to downsampling before transform") |>
  add_option(opt_str = c("-k", "--keycol"),
    type = "character", default = "barcode",
    help = "column of barcode in sa2meta") |>
  add_option(opt_str = c("--dscol"),
    type = "character", default = NULL,
    help = "on which sa2meta column for downsample") |>
  add_option(opt_str = c("--nds"), type = "integer",
    default = 50, help = "# of cells per group for downsample") |>
  add_option(opt_str = c("-n", "--nmax"), type = "integer",
    default = 200000, help = "# cells at most") |>
  parse_args(object = _)

# * DEBUG
datadir <- file.path("/Users/szu/git-recipes/mouseBrainAtlas",
  "amb_pairedtag", "03.integration", "src/test/resource")
args$sa2fnm <- file.path(datadir, "sa2ann_test.h5ad")
args$outdir <- datadir
args$sa2meta <- file.path(datadir, "sa2ann_cellmeta.csv")
args$keycol <- "barcode"
args$dscol <- "brainregion"
args$nds <- 50
args$nmax <- 5000
args$downsample <- FALSE
args$conda <- "~/mambaforge/bin/mamba"
args$condaenv <- "sa2"
args$useBPCells <- TRUE

# * functions
tos <- function(ann,
                outdir,
                meta = NULL) {
  if (!is.null(meta)) {
    message("Cellmeta is not NULL, will save is to seurat too.")
    if (!all(rownames(meta) == colnames(mat))) {
      stop("meta and seurat barcodes are not matched.")
    }
  }
  dir.create(outdir, showWarnings = TRUE)
  outs5file <- file.path(outdir, "anns5.rds")
  message("Transform snapatac2 anndaata to Seurat",
    " with MatrixExtra package.")
  message("Treat X from ann as count, and force it to integer.")
  # ann and seurat rows are different.
  mat <- t(ann$X)
  rownames(mat) <- ann$var_names$to_list()
  colnames(mat) <- ann$obs_names$to_list()
  mat <- MatrixExtra::as.csc.matrix(x = mat) |>
    BPCells::convert_matrix_type(matrix = _, type = "uint32_t")
  ann_meta <- ann$obs |>
    x => `rownames<-`(x, colnames(mat))
  ann_meta$barcode <- rownames(ann_meta)
  if (!is.null(meta)) {
    ann_meta <- merge(ann_meta, meta, by = "barcode") |>
      x => `rownames<-`(x, x$barcode)
  }
  s5 <- Seurat::CreateSeuratObject(counts = mat, meta.data = ann_meta)
  outs5file <- file.path(outdir, "anns5.rds")
  message("Save Seurat to ", outs5file)
  saveRDS(s5, outs5file)
  message("Seurat object saved.")
}

tos5 <- function(ann,
                 outdir,
                 meta = NULL) {
  if (!is.null(meta)) {
    message("Cellmeta is not NULL, will save it to seurat too.")
    if (!all(rownames(meta) == colnames(mat))) {
      stop("meta and seurat barcodes are not matched.")
    }
  }
  dir.create(outdir, showWarnings = TRUE)
  outs5matdir <- file.path(outdir, "_mat")
  dir.create(outs5matdir, showWarnings = TRUE)
  outs5file <- file.path(outdir, "anns5.rds")
  message("Transform snapatac2 anndaata to Seurat5",
    " with BPCells and MatrixExtra package.")
  message("Treat X from ann as count, and force it to integer.")
  # ann and seurat rows are different.
  mat <- t(ann$X)
  rownames(mat) <- ann$var_names$to_list()
  colnames(mat) <- ann$obs_names$to_list()
  message("Save counts to ", outs5matdir, " with BPCells.")
  MatrixExtra::as.csc.matrix(x = mat) |>
    BPCells::convert_matrix_type(matrix = _, type = "uint32_t") |>
    BPCells::write_matrix_dir(mat = _, dir = outs5matdir,
      overwrite = TRUE)
  d <- BPCells::open_matrix_dir(outs5matdir)
  ann_meta <- ann$obs |>
    x => `rownames<-`(x, colnames(mat))
  ann_meta$barcode <- rownames(ann_meta)
  if (!is.null(meta)) {
    ann_meta <- merge(ann_meta, meta, by = "barcode") |>
      x => `rownames<-`(x, x$barcode)
  }
  s5 <- Seurat::CreateSeuratObject(counts = d, meta.data = ann_meta)
  outs5file <- file.path(outdir, "anns5.BPCells.rds")
  message("Save Seurat to ", outs5file)
  saveRDS(s5, outs5file)
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
message("load SnapATAC2 anndata file: ",
  args$sa2fnm)
raw_ann <- sa2$read(filename = args$sa2fnm, backed = "r")

meta <- NULL
if (!is.null(args$sa2meta)) {
  message("load SnapATAC2 meta csv file: ", args$sa2meta,
    " with data.table.")
  meta <- data.table::fread(
    file = args$sa2meta, sep = ",",
    header = TRUE, data.table = FALSE
  ) |> x => `rownames<-`(x, x[[args$keycol]])
  if (args$keycol != "barcode") {
    meta$barcode <- meta[[args$keycol]]
  }
}

if ((!args$downsample) || (raw_ann$n_obs <= args$nmax)) {
  message("Transform ann to seurat without dowmsample")
  # read to memory
  raw_ann <- sa2$read(filename = args$sa2fnm, backed = NULL)
  if (args$useBPCells) {
    tos5(ann = raw_ann, outdir = args$outdir,
      meta = meta, assay = "RNA")
  } else {
    tos(ann = raw_ann, outdir = args$outdir,
      meta = meta, assay = "RNA")
  }
  message("Quit the script. Good Luck!")
  quit(save = "no", status = 0, )
}

if (is.null(args$dscol) || is.null(meta)) {
  message("Downsample before transform.")
  message(str_glue("Downsample to {args$nmax} barcodes."))
  barcodes <- sample(raw_ann$obs_names,
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
  sub_ann <- raw_ann$subset(obs_indices = barcodes,
    out = tmpf)
  sub_ann <- sub_ann$to_memory()
  sub_meta <- meta[barcodes, ]
  if (args$useBPCells) {
    tos5(ann = sub_ann, outdir = args$outdir,
      meta = sub_meta)
  } else {
    tos(ann = sub_ann, outdir = args$outdir,
      meta = sub_meta)
  }
})
message("Quit the script. Good Luck!")
