#' Synthetic data generator from real counts
#'
#' This function generates synthetic count data based on empirical cumulative distribution (ecdf) of real count data
#'
#' @param comm community; a matrix of real count data that we want to simulate/sythesize. Samples are in rows and OTUs are in columns.
#' @param mar MARGIN for apply function to calculate zero proportion for each row (mar = 1) or column (mar = 2).
#' @param Sigma covariance structure of size p by p. p should match with the number of OTUs in \code{comm}, in other words, the number of columns of \code{comm}.
#' @param n number of samples
#' @param seed seed number for data generation (rmvnorm)
#' @param verbose logical value. If it is TRUE, it will print out which iteration is going on and how long it took for calculation for each step. The defulat is FALSE.
#'
#' @return \code{synthData_from_ecdf} returns a data matrix of size n by p.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats ecdf
#' @importFrom stats pnorm
#' @importFrom SpiecEasi make_graph graph2prec prec2cov
#' @export
#'
#' @example man/examples/synthData_ex.R
synthData_from_ecdf <- function(comm, mar = 2, Sigma, n, seed = 10010, verbose = FALSE){

  d <- ncol(comm)
  zratio <- apply(comm, MARGIN = mar, function(x) (sum(x==0)/length(x)))
  maxabund <- apply(comm, MARGIN = mar, max) # to restrict the search range of the solution.

  if(!is.null(seed)) {
    set.seed(seed)
  }
  normd <- mvtnorm::rmvnorm(n, mean=rep(0, d), sigma = Sigma) # mvtnorm package is the fastest one to generate multivariate normal.
  unif <- pnorm(normd)
  dat <- matrix(0, n, d)

  for ( j in 1:d ){
    nzind <- which(unif[, j] > zratio[j]) # to keep the zero ratio as the true data.
    empf <- ecdf(comm[, j]) # empf is a cdf function. empf(c) = Pr(X <= c)

    ptm <- proc.time()
    for ( k in seq_along(nzind) ){
      # This is called "inverse transform sampling". https://en.wikipedia.org/wiki/Inverse_transform_sampling
      # Since the range of the cdf/quantile function is in [0, 1]
      # we want to what is the data value corresponding to a probability (value from unif variable) between 0 and 1.
      # Basically, we numerically calculate the inverse of empirical cdf. find a solution "empf^{-1}(prob)=?"
      dat[nzind[k], j] <- qstepcdf(unif[nzind[k], j], empf, interval = c(0, maxabund[j]))
    }
    if(isTRUE(verbose)) {
      message("iteration = ", j , ": time = ", proc.time() - ptm)
    }

  }
  return(dat)

}

#' Inverse function of empirical cumulative distribution function (ecdf)
#'
#' @param p probability (between 0 and 1)
#' @param empf empirical cdf or any quantile function.
#' @param interval find a solution only within this interval. a vector containing two end points of the interval.
#' @param tol the desired accuracy (convergence tolerance).
#' @param maxiter the maximum number of iterations for \code{uniroot.all}.
#'
#' @importFrom rootSolve uniroot.all
#' @examples
#' ### This is an internal function used in synthData_from_ecdf.
qstepcdf <- function(p, empf, interval, tol = 1e-3, maxiter = 100){
  ans <- c()
  # uniroot.all from rootSolve package was the fastest one.
  sol <- as.numeric(uniroot.all(function(x){empf(x)-p}, interval = interval, tol = tol, maxiter = maxiter))

  # in case quantile function is a step function, added the following step.
  if (p <= empf(floor(sol))) {
    ans <- floor(sol)
  } else if (p > empf(floor(sol))) {
    ans <- ceiling(sol)
  }
  return(ans)
}


