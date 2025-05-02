#!/usr/bin/env Rscript

options(future.globals.maxSize = 10e9)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(optparse)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BSgenome)
library(future.apply)

args <- OptionParser("Iteratively merging peaks from peak summit.
   The method is from ArchR.
") |>
  add_option(
    opt_str = c("-i", "--input"),
    type = "character",
    help = "each line: GroupName\tPeakSummitBedName"
  ) |>
  add_option(
    opt_str = c("-g", "--genome"),
    type = "character", help = "genoem, default is mm10",
    default = "mm10"
  ) |>
  add_option(c("-e", "--extend"),
    type = "integer", default = 500,
    help = "extend bp length, default 500"
  ) |>
  add_option(c("--blacklist"),
    type = "character", default = NULL,
    help = "blacklist bed file, default NULL"
  ) |>
  add_option(c("--chromSize"),
    type = "character", default = NULL,
    help = "chrom size file, default NULL"
  ) |>
  add_option(c("-o", "--out"),
    type = "character", help = "out union peak bed fnm."
  ) |>
  add_option(c("--ncpu"),
    type = "integer",
    default = 4, help = "ncpu for parallel"
  ) |>
  add_option(c("--normalize"), action = "store_true",
    default = FALSE,
    type = "logical",
    help = "if need to noramlize the score") |>
  parse_args(object = _)

plan(multicore, workers = args$ncpu)

# read bed to gr
read2gr <- function(bedF, label) {
  df <- fread(bedF, header = F)
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(
    df$chr,
    IRanges(df$start, df$end)
  )
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}


# extend summit to 500 bp
extendSummit <- function(gr, size = ext_size) {
  gr <- resize(gr, width = size, fix = "center")
  return(gr)
}


# filter blacklist
filter4blacklist <- function(gr, blacklistF) {
  black_list <- read.table(blacklistF)
  black_list.gr <- GRanges(
    black_list[, 1],
    IRanges(black_list[, 2], black_list[, 3])
  )

  idx <- queryHits(
    findOverlaps(gr, black_list.gr)
  )
  if (length(idx) > 0) {
    gr <- gr[-idx]
  }
  return(gr)
}


# filter non-chromosome
filter4chrom <- function(gr, chromF) {
  chrom <- read.table(chromF)
  chrom.gr <- GRanges(
    chrom[, 1],
    IRanges(0, chrom[, 2])
  )
  idx <- queryHits(
    findOverlaps(gr, chrom.gr, type = "within")
  )
  if (length(idx) > 0) {
    gr <- gr[idx]
  }
  return(gr)
}


# filter N containing regions
filter4N <- function(gr, genome = "mm10") {
  genome <- getBSgenome(genome)
  nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))
  mcols(gr)$GC <- round(
    rowSums(nucFreq[, c("G", "C")]) / rowSums(nucFreq), 4
  )
  mcols(gr)$N <- round(nucFreq[, c("N")] / rowSums(nucFreq), 4)
  gr[which(mcols(gr)$N < 0.001)] # Remove N Containing Peaks
  return(gr)
}

#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping
#' set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be
#' used to determine how overlapping regions should be resolved.
#' The resolution of overlapping regions also depends on `decreasing`.
#' For example, if a column named "score" is used for `by`,
#' `decreasing = TRUE` means that the highest "score" in the overlap will be
#' retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the
#' values in the column indicated via `by` should be ordered in decreasing
#' order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value indicating whether the output
#' should include extra reporting.
#' @export
nonOverlappingGR <- function(
    gr = NULL,
    by = "score",
    decreasing = TRUE,
    verbose = FALSE) {
  stopifnot(by %in% colnames(mcols(gr)))
  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE,
                              by = "score", decreasing = TRUE) {
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
    o <- findOverlaps(gr, r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[, by], decreasing = decreasing), ]
    gr <- gr[!duplicated(mcols(gr)$cluster), ]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }

  if (verbose) {
    message("Converging", appendLF = FALSE)
  }
  i <- 0
  grConverge <- gr
  while (length(grConverge) > 0) {
    if (verbose) {
      message(".", appendLF = FALSE)
    }
    i <- i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge,
      filter = TRUE,
      by = by,
      decreasing = decreasing
    )

    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect,
      invert = TRUE,
      ignore.strand = TRUE
    ) # blacklist selected gr

    if (i == 1) { # if i=1 then set gr_all to clustered
      grAll <- grSelect
    } else {
      grAll <- c(grAll, grSelect)
    }
  }
  message(sprintf("Converged after %s iterations!", i))

  if (verbose) {
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))
  return(grAll)
}


# normlize to score per million
norm2spm <- function(gr, by = "score") {
  mlogp <- mcols(gr)[, by]
  normmlogp <- 10^6 * mlogp / sum(mlogp)
  mcols(gr)$spm <- normmlogp
  return(gr)
}


# load summit
summitF <- read.table(args$input, header = F)
label.lst <- as.character(summitF$V1)
file.lst <- as.character(summitF$V2)

peak.list <- future_lapply(
  seq_along(file.lst), function(i) {
    message("working on summit set for... ", label.lst[i])
    p.gr <- read2gr(file.lst[i], label = label.lst[i])
    if (args$extend < 1) {
      message("extend size is less than 1")
      message("No extendsummit is perfomed.")
    } else {
      p.gr <- extendSummit(p.gr, size = args$extend)
    }
    if (!is.null(args$blacklist)) {
      message("filter by blacklist: ", args$blacklist)
      p.gr <- filter4blacklist(p.gr, args$blacklist)
    }
    if (!is.null(args$chromSize)) {
      message("filter by chromSize: ", args$chromSize)
      p.gr <- filter4chrom(p.gr, args$chromSize)
    }
    # p.gr <- filter4N(p.gr, genome = args$genome)
    # p.gr <- nonOverlappingGR(p.gr, by = "score", decreasing = TRUE)
    if(args$normalize) {
      message("perform normalization of score.")
      p.gr <- norm2spm(p.gr, by = "score")
    } else {
      message("treat input score as spm directly.")
      mcols(p.gr)$spm <- mcols(p.gr)[ , "score"]
    }
    return(p.gr)
  }
)
message("Finish filtered peak for each cluster.")

# merge
message("Merge all the peaks.")
merged.gr <- do.call(c, peak.list)
merged.gr <- nonOverlappingGR(merged.gr, by = "spm", decreasing = TRUE)

outUnion <- as.data.frame(merged.gr)

dir.create(dirname(args$out), recursive = TRUE,
  showWarnings = FALSE)
message("save union peak set to: ", args$out)
fwrite(outUnion,
  file = args$out,
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
)
message("Done.")
