#' CMS classification
#' @description consensus molecular subtypes (CMS) classification for mouse tissues.
#' @export
#' @param emat a numeric expression matrix with sample columns, and HGNC symbol rownames.
#' Data should be normalized.
#' see the example in \code{\link{TestData_gemm}}.
#' @param Genesets Select the same as template
#' @param templates a data frame with two columns; \emph{class} (coerced to
#' factor) and \emph{probe} (coerced to character). There are 3 templates: 'template.CMS.A', 'template.CMS.B' and 'template.CMS.C'.
#' Default is template.CMS.A
#' @param nPerm an integer, number of permutations for \eqn{p}-value
#' estimation.
#' @param seed an integer, for \eqn{p}-value reproducibility.
#' @param FDR a false discovery rate, sets prediction confidence threshold.
#' @param verbose a logical, whether console messages are to be displayed.
#' @param doPlot a logical, whether to produce prediction subHeatmap.
#' @details \code{MmCMS} provides 3 options to call CMS subtypes in mouse tissues.
#' The core algorithm is the \href{https://github.com/peterawe/CMScaller}{CMScaller} developed by Eide PW, et al. (2017)
#' and the Nearest Template Prediction (NTP) algorithm as proposed by Yujin Hoshida (2010).
#' @note genes with missing values are discarded.
#' @return a data frame with class predictions, template distances,
#' \eqn{p}-values and false discovery rate adjusted \eqn{p}-values
#' (\code{\link[stats]{p.adjust}}). Rownames equal emat
#' colnames.
#' @seealso \code{\link{template.CMS.A}}, \code{\link{template.CMS.B}}, \code{\link{template.CMS.C}}
#' @references Eide PW, Bruun J, Lothe RA, Sveen A. (2017). CMScaller:
#' an R package for consensus molecular subtyping of colorectal cancer pre-clinical models.
#' doi: 10.1038/s41598-017-16747-x.
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A
#' Single-Sample-Based Flexible Class Prediction with Confidence Assessment.
#' PLoS ONE 5, e15543.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A,
#' Soneson C, et al. The consensus molecular subtypes of colorectal cancer.
#' Nat Med. 2015;21:1350-6.
#' @examples
#' emat <- TestData_gemm
#' re <- MmCMS(emat, templates=MmCMS::template.CMS.A, Genesets = c("template.CMS.A"), seed=120)

MmCMS <- function(emat, templates=MmCMS::template.CMS.A, Genesets = c("template.CMS.B", "template.CMS.C", "template.CMS.A"),
                  nPerm=1000, seed=NULL,
                  FDR=0.05, doPlot=TRUE, verbose=TRUE) {

  # checkInput ##############################################################

  # check datatype input and try to coerce to matrix
  if (inherits(emat, "data.frame")) emat <- as.matrix(emat)
  if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
  if (is.null(rownames(emat))) stop("missing Ensembl id rownames(emat)")

  if (ncol(emat) < 30) warnings("few samples - high prediction variance", call.=FALSE)

  ### Get the GO-BP gene set, saved from MSigDB OR optionC geneset
  if (Genesets == "template.CMS.B") {
    geneset <- GO.BP.list_mus
  } else if (Genesets == "template.CMS.C") {
    geneset <- option.C.geneset.list_mus
  } else {
    geneset <- NULL
  }

  if (!is.null(geneset)) {
    ### Generate single sample GSEA scores (ssGSEA)
    message('Calculating ssGSEA scores...')

    ### set seed
    if (as.Date(paste0(R.version$year,'-',R.version$month,'-',R.version$day)) >= as.Date('2019-04-26')) {
      suppressWarnings(set.seed(123, sample.kind = 'Rounding'))
    } else {
      set.seed(123)
    }

    ### package version change of GSVA
    if (utils::packageVersion("GSVA") <= "1.48.3") {
      ### ssGSEA
      y <- suppressWarnings(GSVA::gsva(emat, geneset, max.sz = Inf, ## min.size default as it may lead to no gene set scores
                                       verbose = F, method = "ssgsea", parallel.sz = 4,
                                       ssgsea.norm = T))
    } else {
      ### ssGSEA -- update
      ssgsea_par <- GSVA::ssgseaParam(emat, geneset,
                                      maxSize = Inf, ## min.size default as it may lead to no gene set scores
                                      normalize = T)
      y <- suppressWarnings(GSVA::gsva(ssgsea_par, verbose = F, BPPARAM = BiocParallel::SnowParam(workers = 2)))
    }

    # scale and center data, basically a wrapper for scale() function
    emat <- ematAdjust(y)
  } else {
    # scale and center data, basically a wrapper for scale() function
    emat <- ematAdjust(emat)
  }

  # ntpPredict ##############################################################
  res <- ntp(emat, templates, seed=seed, nPerm=nPerm, doPlot=doPlot, verbose=verbose)
  res <- subSetNA(res, FDR=FDR, verbose=verbose)

  # output ##################################################################
  # sanity check III - whether any FDR-values are above .1
  if (nPerm > 500 && min(res$FDR) > .1) {
    warning("low-confidence predictions - check input", call.=FALSE)
  }

  return(res)
}
