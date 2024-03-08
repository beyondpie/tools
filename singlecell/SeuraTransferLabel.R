#!/usr/bin/env Rscript

import::from(.from = optparse, make_option,
  OptionParser, add_option, parse_args)
import::from(.from = stringr, str_glue)
## library(Seurat)
## library(tidyverse)
## library(Matrix)
## options(future.globals.maxSize = 5e9)

# This script is to provide a unified transfer label module
# - Only works for two dataset at one time
# - Seurat CCA
# - support co-embedding plot
# - automatically label group based on votes
# - support consensus matrix plot
# - input format is Seurat

# NOTE:
# - R version >= 4.2
# - Seurat version >= 5.0

# * set up arguments
option <- OptionParser() |>
  add_option(opt_str = c("-r", "--ref"), type = "character",
    help = "reference seurat object fnm") |>
  add_option(opt_str = c("-q", "--query"), type = "character",
    help = "query seurat object fnm") |>
  add_option(opt_str = c("-o", "--outdir"), type = "character",
    help = "the directory to save the results") |>
  add_option(opt_str = c("--downsample"), type = "logical",
    action = "store_true", default = FALSE) |>
  add_option(opt_str = c("--drcol"), type = "character",
    help = "col of ref for downsample") |>
  add_option(opt_str = c("--dqcol"), type = "character",
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
  add_option(opt_str = c("-m", "--method"), type = "character",
    default = "cca",
    help = "method for transfer label") |>
  add_option(opt_str = c("-k", "--kanchor"), type = "integer",
    help = "kanchor for finding anchors") |>
  add_option(opt_str = c("-p", "--npca"), type = "integer",
    default = 50, help = "# of components for finding anchors") |>
  add_option(opt_str = c("-f", "--feature"), type = "character",
    default = NULL, help = "gene list file as feature") |>
  add_option(opt_str = c("--nvar"), type = "integer",
    default = 3000, help = " # of variable features") |>
  add_option(opt_str = c("--useref"), type = "logical",
    action = "store_sture", default = TRUE,
    help = "use variable features from ref") |>
  add_option(opt_str = c("--usequery"), type = "logical",
    action = "store_false", default = FALSE) |>
  parse_args(object = _)

# * set up inputs
message("Run Seurat Transfer Label.")
ref <- readRDS(option$ref)
query <- readRDS(option$query)

# * outdir prepare
dir.create(option$outdir, showWarnings = FALSE)
message("Results will be put under: ", option$outdir)

# * downsampling

if (option$downsample) {
  message("Perform downsampling.")
}


# * set up anchor process
# * set up transfer labels and annotations
# * plot co-embedding
# * generate consensus matrix
# * plot consensus matrix
