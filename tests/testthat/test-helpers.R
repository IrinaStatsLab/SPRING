# --- mclr -------------------------------------------------------------------

test_that("mclr returns same dimensions as input", {
  dat <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 2)
  out <- mclr(dat)
  expect_equal(dim(out), dim(dat))
})

test_that("mclr preserves zeros", {
  dat <- matrix(c(0, 1, 0, 2, 3, 0), nrow = 2)
  out <- mclr(dat)
  expect_true(all(out[dat == 0] == 0))
})

test_that("mclr default: all non-zero outputs are strictly positive", {
  dat <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 2)
  out <- mclr(dat)
  expect_true(all(out[dat > 0] > 0))
})

test_that("mclr eps = 0 returns unshifted clr (may have negatives)", {
  dat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  out <- mclr(dat, eps = 0)
  # row means of non-zero log ratios should be zero
  nzero <- dat > 0
  log_dat <- ifelse(nzero, log(dat), 0)
  row_mean_log <- rowMeans(log_dat) / rowMeans(nzero)
  clr <- ifelse(nzero, log_dat - row_mean_log, 0)
  expect_equal(out, clr)
})

test_that("mclr eps > 0 shifts non-zero values by eps", {
  dat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  out0 <- mclr(dat, eps = 0)
  out_shifted <- mclr(dat, eps = 2)
  nzero <- dat > 0
  expect_equal(out_shifted[nzero], out0[nzero] + 2)
  expect_equal(out_shifted[!nzero], rep(0, sum(!nzero)))
})

test_that("mclr works with a single row", {
  dat <- matrix(c(0, 1, 2, 3), nrow = 1)
  out <- mclr(dat)
  expect_equal(dim(out), c(1L, 4L))
  expect_true(all(out[dat == 0] == 0))
  expect_true(all(out[dat > 0] > 0))
})

test_that("mclr accepts a data frame and returns a matrix", {
  dat <- as.data.frame(matrix(c(1, 2, 3, 4), nrow = 2))
  out <- mclr(dat)
  expect_true(is.matrix(out))
})

test_that("mclr warns when atleast is negative and uses default", {
  dat <- matrix(c(1, 2, 3, 4), nrow = 2)
  expect_warning(mclr(dat, atleast = -1), "atleast should be positive")
})

test_that("mclr stops on negative eps", {
  dat <- matrix(c(1, 2, 3, 4), nrow = 2)
  expect_error(mclr(dat, eps = -1))
})

# --- hugeKmb ----------------------------------------------------------------
# Use a small matrix (20 x 5) with 3 lambda values to keep tests fast.

test_that("hugeKmb returns a list with expected elements", {
  data(QMP)
  dat <- QMP[1:20, 1:5]
  lambda <- c(0.5, 0.3, 0.1)
  out <- hugeKmb(dat, lambda = lambda, verbose = FALSE)
  expect_type(out, "list")
  expect_true(all(c("beta", "path", "df", "sparsity", "lambda") %in% names(out)))
})

test_that("hugeKmb beta and path lists have length equal to nlambda", {
  data(QMP)
  dat <- QMP[1:20, 1:5]
  lambda <- c(0.5, 0.3, 0.1)
  out <- hugeKmb(dat, lambda = lambda, verbose = FALSE)
  expect_equal(length(out$beta), length(lambda))
  expect_equal(length(out$path), length(lambda))
})

test_that("hugeKmb beta elements are p x p matrices", {
  data(QMP)
  dat <- QMP[1:20, 1:5]
  p <- ncol(dat)
  lambda <- c(0.4, 0.2)
  out <- hugeKmb(dat, lambda = lambda, verbose = FALSE)
  for (k in seq_along(lambda)) {
    expect_equal(dim(as.matrix(out$beta[[k]])), c(p, p))
  }
})

test_that("hugeKmb sparsity increases with lambda", {
  data(QMP)
  dat <- QMP[1:20, 1:5]
  lambda <- c(0.1, 0.4)
  out <- hugeKmb(dat, lambda = lambda, verbose = FALSE)
  # higher lambda -> sparser graph
  expect_true(out$sparsity[1] >= out$sparsity[2])
})

test_that("hugeKmb works with Rmethod = 'original'", {
  data(QMP)
  dat <- QMP[1:20, 1:5]
  lambda <- c(0.4, 0.2)
  expect_no_error(hugeKmb(dat, lambda = lambda, Rmethod = "original", verbose = FALSE))
})
