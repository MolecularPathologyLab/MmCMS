#' matrix row-wise scaling and centering
#' @export
#' @description Centers and scales gene expression matrix so that each row
#' has mean=0 and sd=1. If input data is non-normalized sequencing \emph{count}
#' data, normalization should be performed by setting \code{normMethod}. This
#' also enforces log2-transformation.
#' @param emat a numeric matrix with row features and sample columns.
#' @param center numeric of same length as \code{nrow(emat)}.
#' @param scale numeric of same length as \code{nrow(emat)}.
#' @param signalFilt numeric setting feature filtering. Specifically rows
#' where \code{rowMax(emat) < quantile(emat, probs = signalFilt)} are discarded.
#' @param verbose logical, console messages output.
#' @param ... additional arguments passed to normalizeBetweenArrays and
#' calcNormFactors depending normalization method selected.
#' @return
#' A row-wise centered and scaled matrix. Output matrix will have fewer rows
#' than input if signalFilt>0.
#' @details
#' \code{ematAdjust} performs row-wise scaling and centering by passing matrix
#' to \code{\link[base]{scale}}. Setting \code{scale} and \code{center} may be
#' useful for predicting new samples based on row-wise means and standard
#' deviations from prior (identically processed) datasets. If \emph{e.g.}
#' \code{signalFilt=.1} features with maximum below emat 10th
#' percentile are discarded.
ematAdjust <- function(emat, center = TRUE, scale = TRUE,
                       signalFilt = 0, verbose = getOption("verbose"), ...)
{

  # checkInput ##############################################################
  if (is.data.frame(emat)) emat <- as.matrix(emat)

  N <- ncol(emat)
  P.in <- nrow(emat)

  # processData #############################################################

  # filter low signal probes
  if (signalFilt>0) {
    signal.filter <- stats::quantile(emat, signalFilt, na.rm = TRUE)
    filterLow <- apply(emat,1, max, na.rm = TRUE) < signal.filter
    emat <- emat[!filterLow,]
    if (length(center) > 1) center <- center[!filterLow]
    if (length(scale) > 1) center <- scale[!filterLow]
  } else {
    filterLow=FALSE
  }

  # standardize
  emat <- t(scale(t(emat), scale=scale, center=center))

  P.out <- nrow(emat)
  isnorm <- NULL

  if (verbose == TRUE) {
    emat.sd <- round(stats::sd(emat, na.rm = TRUE),2)
    emat.mean <- round(mean(emat,na.rm = TRUE),2)
    if (abs(emat.mean) >.5) { isnorm = " <- check processing!" }
    message(paste0("input: ", N, " samples; ", P.in, " features"))
    message(paste0("output: ", N, " samples; ", P.out, " features; ",
                   "mean=", emat.mean, "; sd=", emat.sd, isnorm))}

  # returnData ##############################################################
  return(emat)
}
