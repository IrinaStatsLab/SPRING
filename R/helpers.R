# Rfunctions for simple example for SPRING method
# Yoon, Gaynanova and M\"{u}eller (2019) Frontiers in Genetics, Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data.
# doi:10.3389/fgene.2019.00516

# To implement Kendall correlation estimates on huge.
# original huge function can only take covariance or data matrix.
# huge.mb funciton returns beta values (MB coefficient estimates)
# huge function with method="mb" does not return beta values.
# to do network visualization (for edge color, I need beta)


#' Internal wrapper function to implement rank-based correlation to huge.mb function in "huge" package
#'
#' @param data n by p matrix data. usually through pulsar, data will receive subsamples.
#' @param lambda a vector of lambda values
#' @param type a type of variables passed to \code{latentcor}. "tru" (truncated) is default. See \code{latentcor::latentcor} for valid values ("con", "bin", "tru", "ter").
#' @param sym "or" is the symmetrizing rule of the output graphs. If sym = "and", the edge between node i and node j is selected ONLY when both node i and node j are selected as neighbors for each other. If sym = "or", the edge is selected when either node i or node j is selected as the neighbor for each other. The default value is "or". (refer to huge manual)
#' @param verbose If \code{verbose = FALSE}, tracing information printing for HUGE (High-dimensional Undirected Graph Estimation) is disabled. The default value is TRUE.
#' @param Rmethod The calculation method of latent correlation. Either "approx" or "original". If \code{Rmethod = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{Rmethod = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD Logical indicator. \code{use.nearPD = TRUE} gets nearest positive definite matrix for the estimated latent correlation matrix with shrinkage adjustment by \code{nu}. Output \code{R} is the same as \code{Rpointwise} if \code{use.nearPD = FALSE}. Default value is \code{TRUE}.
#' @param nu Shrinkage parameter for the correlation matrix, must be between 0 and 1. Guarantees that the minimal eigenvalue of the returned correlation matrix is greater or equal to \code{nu}. The default (recommended) value is 0.001.
#' @param ratio When \code{Rmethod = "approx"}, specifies the boundary value for multilinear interpolation, must be between 0 and 1. The default (recommended) value is 0.9. Ignored when \code{Rmethod = "original"}.
#' @param tol When \code{Rmethod = "original"}, desired accuracy of the bridge function inversion. Ignored when \code{Rmethod = "approx"}. Default is 1e-6.
#'
#' @return \code{hugeKmb} returns a list containing
#' \describe{
#'      \item{beta}{a list of length \code{nlambda}, each element is a sparse p by p matrix of MB coefficient estimates at the corresponding lambda value.}
#'      \item{path}{a list of length \code{nlambda}, each element is a sparse p by p adjacency matrix of selected edges (symmetrized) at the corresponding lambda value.}
#'      \item{df}{a vector of length \code{nlambda * p} giving the number of selected neighbors for each node at each lambda value.}
#'      \item{sparsity}{a vector of length \code{nlambda} giving the proportion of selected edges at each lambda value.}
#'      \item{lambda}{the lambda sequence used.}
#'      \item{cov.input}{logical, whether the input was identified as a covariance matrix.}
#'      \item{scr}{logical, whether screening was applied.}
#' }
#'
#' @importFrom huge huge.mb
#' @export
#'
#' @examples
#' data(QMP)
#' lambda <- c(0.5, 0.3, 0.1)
#' out <- hugeKmb(QMP[, 1:10], lambda = lambda, verbose = FALSE)
#'
hugeKmb <- function(data, lambda, type = "tru", sym = "or", verbose = TRUE, Rmethod = c("approx", "original"), tol = 1e-6, use.nearPD = TRUE, nu = 0.001, ratio = 0.9) {
  Rmethod <- match.arg(Rmethod)
  S    <- latentcor::latentcor(data, types = type, method = Rmethod, tol = tol, use.nearPD = use.nearPD, nu = nu, ratio = ratio, showplot = FALSE)$R
  est  <- huge::huge.mb(S, lambda, sym = sym, verbose = verbose)
  est
}




#' Modified central log ratio (mclr) transformation
#'
#' @param dat raw count data or compositional data (n by p) does not matter.
#' @param base exp(1) for natural log
#' @param tol tolerance for checking zeros

# For eps and atleast, users do not have to specify any values. Default should be enough.
#' @param eps epsilon in eq (2) of the paper "Yoon, Gaynanova, Müller (2019), Frontiers in Genetics". positive shifts to all non-zero compositions. Refer to the paper for more details. eps = absolute value of minimum of log ratio counts plus c.
#' @param atleast default value is 1. Constant c which ensures all nonzero values to be strictly positive. default is 1.
#'
#'
#' @return \code{mclr} returns a data matrix of the same dimension with input data matrix.
#' @export
#'
#' @examples
#' data(QMP)
#' RMP <- QMP/rowSums(QMP)
#' mclr_RMP <- mclr(RMP)
#'
mclr <- function(dat, base = exp(1), tol = 1e-16, eps = NULL, atleast = 1){
  dat <- as.matrix(dat)
  nzero <- (dat >= tol)  # index for nonzero part
  LOG <- ifelse(nzero, log(dat, base), 0.0) # take log for only nonzero values. zeros stay as zeros.

  # centralize by the log of "geometric mean of only nonzero part" # it should be calculated by each row.
  if (nrow(dat) > 1){
    clrdat <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0.0)
  } else if (nrow(dat) == 1){
    clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
  }

  if (is.null(eps)){
    if(atleast < 0){
      warning("atleast should be positive. The functions uses default value 1 instead.")
      atleast <- 1
    }
    if( min(clrdat) < 0 ){ # to find the smallest negative value and add 1 to shift all data larger than zero.
      positivecst <- abs(min(clrdat)) + atleast # "atleast" has default 1.
    }else{
      positivecst <- 0
    }
    # positive shift
    ADDpos <- ifelse(nzero, clrdat + positivecst, 0.0) ## make all non-zero values strictly positive.
    return(ADDpos)
  } else if(eps == 0){
    ## no shift. clr transform applied to non-zero proportions only. without pseudo count.
    return(clrdat)
  } else if(eps > 0){
    ## use user-defined eps for additional positive shift.
    ADDpos <- ifelse(nzero, clrdat + eps, 0.0)
    return(ADDpos)
  } else {
    stop("check your eps value for additional positive shift. Otherwise, leave it as NULL.")
  }
}
