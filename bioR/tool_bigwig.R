#' Get bigwig signals for a given genomic ranges
#'
#' @param bwGR use rtracklayer::import to load bigwig file directly.
#' The results will be genomic ranges object.
#' @param peakGR genomic ranges interested
#' We assumme a name field for the ranges, which will be the names
#' for the output vector.
#' @param fnbw how to handle the overlapped bigwig signal itself
#' fn: Double => Double, default identity
#' @param fnGroup how to handle all the overlapped bigwig signals for
#' a given region. fn: Vector[Double] => Double, default mean
#' @return named numeric vector
#' The order the vector follows the input of peakGR
#' @export
getbwSignal <- function(bwGR, peakGR, fnbw = identity, fnGroup = mean) {
  ovlp <- GenomicRanges::findOverlaps(query = peakGR, subject = bwGR)
  s1 <- S4Vectors::mcols(bwGR[S4Vectors::subjectHits(ovlp)])$score |>
    x => fnbw(x)
  s2 <- S4Vectors::aggregate(
    s1,
    by = list(S4Vectors::queryHits(ovlp)), FUN = fnGroup
  )
  r <- rep(NA, length(peakGR))
  names(r) <- S4Vectors::mcols(peakGR)$name
  peaknms <- S4Vectors::mcols(peakGR[s2[, 1]])$name
  r[S4Vectors::mcols(peakGR[s2[, 1]])$name] <- s2[, 2]
  return(r)
}

#' Load Bigwig file into a Genomic range object.
#' @param bwfnm bigwig file name
#' @export
loadBigWigGR <- function(bwfnm) {
  rtracklayer::import(con = bwfnm)
}
