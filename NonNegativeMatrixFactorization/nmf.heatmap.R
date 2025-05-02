# Draw Heatmap of NMF.
library(optparse)
library(ComplexHeatmap)

op <- list(
  make_option(c("--cpmCbyPh5"), type = "character"),
  make_option(c("--peakStatWfnm"), type = "character"),
  make_option(c("--figd"), type = "character"),
  make_option(c("--module"), type = "integer", default = 150),
  make_option(c("--tag"), type = "character", default = "ppdc")
)
args <- parse_args(OptionParser(option_list = op))

# debug
## args$cpmCbyPh5 <- file.path(
##   "/projects/ps-renlab2/szu/projects/amb_pairedtag/",
##   "data/ptHistoneCounts/ATACPeak_rds",
##   "cpm.scbyp.H3K27ac.h5"
## )
## args$peakStatWfnm <- file.path(
##   "/projects/ps-renlab2/szu/projects/amb_pairedtag/",
##   "05.CRE/out/nmf",
##   "nmfPmat.H3K27ac.r150.n0.statW"
## )
## args$figd <- file.path(
##   "/projects/ps-renlab2/szu/projects/amb_pairedtag",
##   "05.CRE/out/nmf"
## )
## args$module <- 150
## args$tag <- "H3K27ac"


# * configs
module <- args$module
tag <- args$tag
statPeakFile <- args$peakStatWfnm
peakBed <- args$peakBed
clusterOrd <- args$clusterOrd
cpmPeakByCluster <- args$cpmCbyPh5

# * functions
loadStatPeakNMF <- function(fnm) {
  statPeak <- data.table::fread(
    file = fnm,
    header = FALSE, sep = "\t", data.table = FALSE
  )
  colnames(statPeak) <- c(
    "peak", "index", "class0", "contributes",
    "featureScore", "selt_fs_list",
    "selt_med_list")
  statPeak$moduleN <- paste(
    "m", as.integer(statPeak$class0 + 1), sep = "")
  rownames(statPeak) <- statPeak$peak
  return(statPeak)
}
getPeakModuleAvgScore.NMF <- function(cpm, moduleOfPeak,
                                      fn = rowMeans) {
  modules <- unique(moduleOfPeak)
  avgScores <- vapply(modules, function(m) {
    a <- cpm[, moduleOfPeak %in% m, drop = FALSE]
    return(fn(a))
  }, FUN.VALUE = rep(0, nrow(cpm)))
  return(avgScores)
}

getTopRankClusterForPeakModule <- function(cluster2moduleScore,
                                           modules = colnames(cluster2moduleScore),
                                           topk = 3) {
  n_module <- ncol(cluster2moduleScore)
  r <- vapply(seq_len(n_module), function(i){
    return(
      order(cluster2moduleScore[, i], decreasing = TRUE)[1:topk]
    )
  },FUN.VALUE = rep(1, topk))
  if(topk < 2) {
    r <- matrix(data = r, nrow = 1)
  }
  colnames(r) <- modules
  return(r)
}

rankModules <- function(topRank.module,
                        modules = colnames(topRank.module),
                        avg.fn = colMeans) {
  avgRank <- avg.fn(topRank.module)
  index.order <- order(avgRank, decreasing = FALSE)
  return(modules[index.order])
}

colFn <- function(col.fn = min) {
  a <- function(mat) {
    apply(mat, 2, col.fn)
  }
  return(a)
}

orderModules.cpm.statPeak <- function(cpm, statPeak,
                                      fn.avgModule = rowMeans,
                                      topk = 1,
                                      fn.rank = min) {
  moduleOfPeak <- statPeak$module
  cluster2moduleScore <- getPeakModuleAvgScore.NMF(
    cpm = cpm,
    moduleOfPeak = moduleOfPeak,
    fn = fn.avgModule)
  topRank.module <- getTopRankClusterForPeakModule(
    cluster2moduleScore = cluster2moduleScore,
    modules = colnames(cluster2moduleScore),
    topk = topk
  )
  colMins <- colFn(col.fn = fn.rank)
  modules.order <- rankModules(
    topRank.module,
    modules = colnames(topRank.module),
    avg.fn = colMins)
  return(modules.order)
}

orderStatPeakByModuleRank <- function(statPeak, mod.ordered){
  index.list <- lapply(mod.ordered, function(m) {
    index.raw <- which(statPeak$moduleN %in% m)
    index <- index.raw[
      order(statPeak$featureScore[index.raw], decreasing = TRUE)]
    return(index)
  })
  index.ordered <- unlist(index.list)
  r <- statPeak[index.ordered, ]
  return(r)
}

downsample.index.1 <- function(index.all,
                               mod.index,
                               mod.ordered = NULL,
                               score.index = NULL,
                               size = 200,
                               seed = 2022) {
  set.seed(seed = seed)
  mods <- if(!is.null(mod.ordered)) {
    mod.ordered
  } else {
    unique(mod.index)
  }
  index.list <- lapply(mods, function(m) {
    index <- mod.index %in% m
    if(sum(index) < 1) {
      warning("Module ", m, " has no features.")
      return(NULL)
    }
    cols <- which(index)
    cols.sampled <- if(!is.null(score.index)) {
      cols[order(score.index[cols], decreasing = TRUE)][
        seq_len(min(length(cols), size))
      ]
    } else {
      sample(cols,size = min(length(cols), size),
        replace = FALSE)
    }
    return(cols.sampled)
  })
  index.list[sapply(index.list, is.null)] <- NULL
  return(unlist(index.list))
 }

downsampleStatPeak.1 <- function(statPeak,
                                 mod.ordered = NULL,
                                 size = 200,
                                 seed = 2022) {
  index.sampled <- downsample.index.1(
    index.all = seq_len(nrow(statPeak)),
    mod.index = statPeak$moduleN,
    mod.ordered = mod.ordered,
    score.index = statPeak$featureScore,
    size = size,
    seed = seed
  )
  statPeak.sampled <- statPeak[index.sampled, ]
  return(statPeak.sampled)
}

getNMFHeatMap <- function(cpm.plot,
                          ha_row = NULL,
                          ha_col = NULL,
                          showRowNames = TRUE,
                          showColNames = FALSE,
                          fontsize = 6,
                          low.val.col = quantile(cpm.plot, 0.01),
                          high.val.col = quantile(cpm.plot, 0.99),
                          use_raster = TRUE,
                          legend_title = "log of CPM",
                          legend_labels = c("Low", "High")) {
  col_fun <- circlize::colorRamp2(
    seq(low.val.col, high.val.col, length = 60),
    viridis::viridis(60)
  )
  p <- ComplexHeatmap::Heatmap(
    matrix = cpm.plot,
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = showRowNames,
    row_names_gp = grid::gpar(fontsize = fontsize),
    show_column_names = showColNames,
    top_annotation = ha_col,
    left_annotation = ha_row,
    use_raster = use_raster,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = legend_title,
      at = c(low.val.col, high.val.col),
      labels = legend_labels)
  )
  return(p)
}

# ============== Main ================
# * load cpm, subclass order and colors
fh5 <- hdf5r::H5File$new(cpmPeakByCluster, mode = "r")
cpm <- fh5[["X/mat"]][,]
peakCoords <- fh5[["X/colnames"]][]
cluster.order.hc <- fh5[["X/rownames"]][]
rownames(cpm) <- cluster.order.hc
colnames(cpm) <- peakCoords

## TODO: update
ha_row <- NULL

# * load nmf result
statPeak <- loadStatPeakNMF(statPeakFile)
n_total <- 10000
n_module <- length(unique(statPeak$moduleN))
n_avg <- floor(n_total / n_module)
statPeak.ds <- downsampleStatPeak.1(
  statPeak = statPeak,
  mod.ordered = NULL,
  size = n_avg,
  seed = 2024)

# scale(x = _, center = TRUE, scale = TRUE) |>
logcpm <- log1p(cpm) |>
  x => x[, statPeak.ds$peak]

# * set module order
nmf.order <- orderModules.cpm.statPeak(
  cpm = logcpm,
  statPeak = statPeak.ds)

# * re-order statPeak
statPeak.ordered <- orderStatPeakByModuleRank(
  statPeak = statPeak.ds,
  mod.ordered = nmf.order)

# * reorder peaks by putting the global modules firstly
modules <- unique(statPeak.ordered$moduleN)
cpm.mod <- getPeakModuleAvgScore.NMF(
  cpm = logcpm,
  moduleOfPeak = statPeak.ds$moduleN,
  fn = rowMeans
)
high.value <- quantile(logcpm, 0.95)
nsubclass.peak.plot <- colSums(cpm.mod >= high.value)
global_modules <- colnames(cpm.mod)[
   which(nsubclass.peak.plot >= (0.7 * nrow(logcpm)))]
mod.reorder <- c(
  global_modules, setdiff(nmf.order, global_modules)
)
statPeak.reorder <- orderStatPeakByModuleRank(
  statPeak.ordered, mod.reorder
)


# * set heatmap column annotation
moduleColor <- grDevices::colorRampPalette(
  colors = c("bisque", "burlywood4") )(
    length(mod.reorder)) |>
    setNames(object = _, mod.reorder)

ha_col <- ComplexHeatmap::HeatmapAnnotation(
  Module = statPeak.reorder$moduleN,
  FeatureScore = ComplexHeatmap::anno_barplot(
    statPeak.reorder$featureScore),
  col = list(
    Module = moduleColor
  ),
  which = "column",
  show_legend = c(FALSE))

# * get heatmap
p.cpm.log2 <- getNMFHeatMap(
  cpm.plot = logcpm[ , statPeak.reorder$peak],
  ## NOTE: update later
  ha_row = NULL,
  ha_col = ha_col,
  fontsize = 6,
  low.val.col = quantile(logcpm, 0.01),
  high.val.col = quantile(logcpm, 0.99))

# * plot figures
figfnm <- file.path(args$figd,
    paste(tag, paste0("r", module), "pdf", sep = "."))
withr::with_pdf(
  new = figfnm,
  code = {
    print(p.cpm.log2)
  }, height = 28, width = 20
)
message("plot figure done: ", figfnm)

# * save meta data
saveRDS(
  object = list(logCPM4plot = logcpm[ , statPeak.reorder$peak],
    nmfModule4plot = statPeak.reorder
  ),
  file = file.path(args$figd,
    paste(tag, paste0("r", module), "meta.rds", sep = "."))
)
