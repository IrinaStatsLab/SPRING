#' Semi-Parametric Rank-based approach for INference in Graphical model (SPRING)
#'
#' @description SPRING follows the neighborhood selection methodology outlined in "mb" method (Meinshausen and Buhlmann (2006)).
#'
#' @param data n by p matrix of microbiome count data, either quantitative or compositional counts. Each row represents each subject/sample and each column represents each OTU (operational taxonomic unit).
#' @param quantitative default is FALSE, which means input "data" is compositional data, which will be normalized using mclr transformation within a function. If TRUE, it means "quantitative" counts are input and no normalization will be applied.
#' @param method graph estimation method. Currently only \code{"mb"} (Meinshausen-Bühlmann neighborhood selection) is available.
#' @param lambda.min.ratio ratio of the smallest to largest value in the lambda sequence. Default is 0.01.
#' @param nlambda number of lambda values in the regularization sequence. Default is 20.
#' @param lambdaseq a sequence of decreasing positive numbers to control the regularization. The default sequence has 20 values generated to be equally spaced on a logarithmic scale from 0.6 to 0.006. Users can specify a sequence to override the default sequence. If user specify as "data-specific", then the lambda sequence will be generated using estimated rank-based correlation matrix from data.
#' @param seed the seed for subsampling.
#' @param ncores number of cores to use for subsampling. The default is 1.
#' @param thresh threshold for StARS selection criterion. 0.1 is recommended (default). The smaller threshold returns sparser graph.
#' @param subsample.ratio 0.8 is default. The recommended values are 10*sqrt(n)/n for n > 144 or 0.8 otherwise.
#' @param rep.num the repetition number of subsampling for StARS edge stability selection. The default value is 20.
#' @param Rtol Desired accuracy when calculating the solution of bridge function in latentcor. Default is 1e-6.
#' @param verbose If \code{verbose = FALSE}, tracing information printing for HUGE (High-dimensional Undirected Graph Estimation) is disabled. The default value is TRUE.
#' @param Rmethod The calculation method of latent correlation. Either \code{"approx"} or \code{"original"}. If \code{Rmethod = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{Rmethod = "original"}, optimization of the bridge inverse function is used. The default is \code{"approx"}.
#' @param use.nearPD Logical indicator. \code{use.nearPD = TRUE} gets nearest positive definite matrix for the estimated latent correlation matrix with shrinkage adjustment by \code{nu}. Output \code{R} is the same as \code{Rpointwise} if \code{use.nearPD = FALSE}. Default value is \code{TRUE}.
#' @param nu Shrinkage parameter for the correlation matrix, must be between 0 and 1. Guarantees that the minimal eigenvalue of the returned correlation matrix is greater or equal to \code{nu}. The default (recommended) value is 0.001.
#' @param ratio When \code{Rmethod = "approx"}, specifies the boundary value for multilinear interpolation, must be between 0 and 1. The default (recommended) value is 0.9. Ignored when \code{Rmethod = "original"}.
#'
#' @return \code{SPRING} returns a list containing
#' \itemize{
#'       \item{output: }{Output results of \code{pulsar::pulsar} based on StARS criterion. It contains:}
#'               \itemize{
#'                        \item{merge: } a list of length \code{nlambda} and each element of list contains a matrix of edge selection probability. Each lambda value, this edge selection probability is calculated across \code{rep.num}.
#'                        \item{summary: } the summary statistic over \code{rep.num} graphs at each value of lambda
#'                        \item{opt.index: } index (along the path) of optimal lambda selected by the criterion at the desired threshold. Will return \eqn{0} if no optimum is found or \code{NULL} if selection for the criterion is not implemented.
#'                        \item{criterion: } we use StARS for our stability criterion.
#'               }
#'       \item{fit: }{Output results of \code{pulsar::refit} function. It contains:}
#'       \itemize{
#'               \item{est: } a list containing
#'               \itemize{
#'                        \item{beta: } Estimates of beta coefficient matrices (of size p by p) by "mb" method on the whole data at each of whole lambda sequence value.
#'                        \item{path: } Estimates of precision matrix (of size p by p) on the whole data at each of whole lambda sequence value.
#'               }
#'               \item{refit: } final estimates of precision matrix (of size p by p).
#'               }
#'       \item{lambdaseq: }{lambda sequence used in the analysis}
#' }
#' @importFrom pulsar pulsar
#' @importFrom latentcor latentcor
#'
#' @export
#'
#' @references
#'
#' Meinshausen N. and Buhlmann P. (2006) \href{https://projecteuclid.org/download/pdfview_1/euclid.aos/1152540754}{"High-dimensional graphs and variable selection with the lasso"}, \emph{The Annals of Statistics}, Vol 34, No. 3, 1436 - 1462.
#'
#' Yoon G., Gaynanova I. and Müller C. (2019) \href{https://www.frontiersin.org/articles/10.3389/fgene.2019.00516/full}{"Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data"}, \emph{Frontiers in Genetics}, 10:516.
#'
#' @example man/examples/ex.R
#'
SPRING <- function(data, quantitative = FALSE, method = c("mb"), lambda.min.ratio = 1e-2, nlambda = 20, lambdaseq = exp(seq(log(0.6), log(0.6*lambda.min.ratio), length.out = nlambda)), seed = 10010, ncores = 1, thresh = 0.1, subsample.ratio = 0.8, rep.num = 20, Rtol = 1e-6, verbose = TRUE, Rmethod = "approx", use.nearPD = TRUE, nu = 0.001, ratio = 0.9){

  method <- match.arg(method)

  if (any(data < 0)) {
    stop("Negative values are detected, but either quantitative or compositional counts are expected.\n")
  }
  p <- ncol(data)
  if (quantitative){
    if (max(rowSums(data)) <= 1 | isTRUE(all.equal(max(rowSums(data)), 1))){
      warning("The input data is normalized, but quantitative count data is expected.\n")
    }
    qdat <- data
  } else {
    qdat <- mclr(data)
  }
  rm(data)
  gc()

  if(is.character(lambdaseq)){
    if(lambdaseq == "data-specific"){
      Kcor <- latentcor::latentcor(qdat, types = "tru", method = Rmethod, tol = Rtol, use.nearPD = use.nearPD, nu = nu, ratio = ratio)$R
      # generate lambda sequence
      lambda.max <- max(max(Kcor-diag(p)), -min(Kcor-diag(p)))
      lambda.min <- lambda.min.ratio * lambda.max
      lambdaseq <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    } else {
      stop("The input for lambdaseq is not correct.\n")
    }
  }

  if(method == "mb"){
    fun <- hugeKmb
  }

  out1.K_count <- pulsar::pulsar(qdat, fun = fun, fargs = list(lambda = lambdaseq, Rmethod = Rmethod, tol = Rtol, verbose = verbose, use.nearPD = use.nearPD, nu = nu, ratio = ratio), rep.num = rep.num, criterion = 'stars', seed = seed, ncores = ncores, thresh = thresh, subsample.ratio = subsample.ratio)

  fit1.K_count <- pulsar::refit(out1.K_count)

  return(list(output = out1.K_count, fit = fit1.K_count, lambdaseq = lambdaseq))
}
