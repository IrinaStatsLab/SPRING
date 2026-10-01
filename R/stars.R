# StARS (Stability Approach to Regularization Selection) for SPRING.
#
# Liu H., Roeder K. and Wasserman L. (2010) "Stability Approach to
# Regularization Selection (StARS) for High Dimensional Graphical Models",
# Advances in Neural Information Processing Systems 23.
#
# This reimplements the subset of pulsar::pulsar() and pulsar::refit() used by
# SPRING (criterion = "stars", no lower/upper bound shortcuts), so results are
# identical to SPRING versions that depended on the pulsar package.

# Run StARS over the lambda path.
#
# data: n by p matrix passed to `fun` (and subsampled by row).
# fun: estimator; called as do.call(fun, c(fargs, list(data))) and must return
#   a list with member `path` (one p by p adjacency matrix per lambda).
# fargs: list of arguments for `fun`; must contain `lambda`.
#
# Returns a list with `stars` (merge, summary, opt.index, criterion, thresh)
# and `est` (the output of `fun` on the full data).
.stars <- function(data, fun, fargs, thresh = 0.1, subsample.ratio = NULL,
                   rep.num = 20, seed = NULL, ncores = 1) {
  n <- nrow(data)
  p <- ncol(data)
  nlams <- length(fargs$lambda)
  if (is.null(subsample.ratio)) {
    subsample.ratio <- if (n > 144) 10 * sqrt(n) / n else 0.8
  }

  # Draw the subsamples up front so results do not depend on `ncores`.
  # The caller's RNG state is restored afterwards.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old.seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old.seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  ind.sample <- replicate(rep.num, sample(seq_len(n), floor(n * subsample.ratio), replace = FALSE), simplify = FALSE)

  estFun <- function(ind) {
    do.call(fun, c(fargs, list(data[ind, , drop = FALSE])))$path
  }

  # Full-data fit, used for the refit at the selected lambda.
  fullmodel <- do.call(fun, c(fargs, list(data)))
  if (!("path" %in% names(fullmodel))) {
    stop("expected data structure with 'path' member")
  }

  paths <- .lapply_cores(ind.sample, estFun, ncores = ncores)
  errors <- vapply(paths, inherits, logical(1), what = "try-error")
  if (all(errors)) {
    stop("all ", rep.num, " subsamples failed with: ", paths[[1]], call. = FALSE)
  }
  if (any(errors)) {
    warning(sum(errors), " of ", rep.num, " subsamples failed and were dropped: ", paths[errors][[1]], call. = FALSE)
    paths <- paths[!errors]
    rep.num <- length(paths)
  }

  # Edge selection frequency across subsamples, for each lambda
  merge <- lapply(seq_len(nlams), function(i) {
    Reduce("+", lapply(paths, function(path) path[[i]])) / rep.num
  })
  # Total instability, made monotone along the (decreasing) lambda path
  summary <- vapply(merge, function(m) 4 * sum(m * (1 - m)) / (p * (p - 1)), numeric(1))
  summary <- cummax(summary)

  if (any(summary >= thresh)) {
    # Last lambda before instability first reaches `thresh`
    opt.index <- max(which.max(summary >= thresh) - 1, 1)
  } else {
    warning("Optimal lambda may be outside the supplied values")
    opt.index <- nlams
  }
  if (opt.index == 1) {
    direction <- if (any(summary >= 0.1)) "larger" else "smaller"
    warning("Optimal lambda may be ", direction, " than the supplied values")
  }

  list(
    stars = list(merge = merge, summary = summary, opt.index = opt.index, criterion = "stars.stability", thresh = thresh),
    est = fullmodel
  )
}

# Full-data estimates and the graph at the StARS-selected lambda
.stars_refit <- function(obj) {
  list(est = obj$est, refit = list(stars = obj$est$path[[obj$stars$opt.index]]))
}

# lapply over subsamples, optionally on a PSOCK cluster (works on all
# platforms). Errors are returned as "try-error" strings, not raised.
.lapply_cores <- function(X, FUN, ncores = 1) {
  safeFUN <- function(x) {
    tryCatch(FUN(x), error = function(e) structure(conditionMessage(e), class = "try-error"))
  }
  if (ncores > 1) {
    cl <- parallel::makePSOCKcluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapplyLB(cl, X, safeFUN)
  } else {
    lapply(X, safeFUN)
  }
}
