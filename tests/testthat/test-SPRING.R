# Helper: run SPRING on a small subset with minimal iterations.
# suppressWarnings silences pulsar's "Optimal lambda may be larger than
# supplied values" which occurs with small p and few lambdas in tests.
spring_small <- function(...) {
  suppressWarnings(
    SPRING(QMP[, 1:10],
      quantitative = TRUE, nlambda = 5, rep.num = 5,
      verbose = FALSE, ...
    )
  )
}

# --- synthData_from_ecdf ----------------------------------------------------
# Use a small subset of QMP (5 taxa) to keep tests fast.

test_that("synthData_from_ecdf returns matrix of correct dimensions", {
  data(QMP)
  comm <- QMP[, 1:5]
  Sigma <- diag(5)
  set.seed(1)
  out <- synthData_from_ecdf(comm, Sigma = Sigma, n = 20)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(20L, 5L))
})

test_that("synthData_from_ecdf returns non-negative values", {
  data(QMP)
  comm <- QMP[, 1:5]
  Sigma <- diag(5)
  set.seed(1)
  out <- synthData_from_ecdf(comm, Sigma = Sigma, n = 20)
  expect_true(all(out >= 0))
})

test_that("synthData_from_ecdf is reproducible with the same seed", {
  data(QMP)
  comm <- QMP[, 1:5]
  Sigma <- diag(5)
  set.seed(42)
  out1 <- synthData_from_ecdf(comm, Sigma = Sigma, n = 10)
  set.seed(42)
  out2 <- synthData_from_ecdf(comm, Sigma = Sigma, n = 10)
  expect_equal(out1, out2)
})

test_that("synthData_from_ecdf preserves zero structure", {
  data(QMP)
  comm <- QMP[, 1:5]
  Sigma <- diag(5)
  set.seed(1)
  out <- synthData_from_ecdf(comm, Sigma = Sigma, n = 50)
  # zero ratio in output should be in the same ballpark as input (within 20%)
  zratio_in <- apply(comm, 2, function(x) mean(x == 0))
  zratio_out <- apply(out, 2, function(x) mean(x == 0))
  expect_true(all(abs(zratio_out - zratio_in) < 0.20))
})

# --- SPRING -----------------------------------------------------------------
# Use a small subset of QMP (10 taxa) with few lambda/subsamples for speed.

test_that("SPRING returns a list with output, fit, lambdaseq", {
  data(QMP)
  fit <- spring_small()
  expect_type(fit, "list")
  expect_true(all(c("output", "fit", "lambdaseq") %in% names(fit)))
})

test_that("SPRING lambdaseq has correct length", {
  data(QMP)
  fit <- spring_small()
  expect_equal(length(fit$lambdaseq), 5)
})

test_that("SPRING opt.index is within valid range", {
  data(QMP)
  fit <- spring_small()
  opt <- fit$output$stars$opt.index
  expect_true(opt >= 1 && opt <= 5)
})

test_that("SPRING adjacency matrix is symmetric", {
  data(QMP)
  fit <- spring_small()
  adj <- as.matrix(fit$fit$refit$stars)
  expect_equal(adj, t(adj))
})

test_that("SPRING works with compositional data (quantitative = FALSE)", {
  data(QMP)
  RMP <- QMP[, 1:10] / rowSums(QMP[, 1:10])
  expect_no_error(suppressWarnings(
    SPRING(RMP, quantitative = FALSE, nlambda = 5, rep.num = 5, verbose = FALSE)
  ))
})

test_that("SPRING works with data-specific lambda sequence", {
  data(QMP)
  fit <- spring_small(lambdaseq = "data-specific")
  expect_equal(length(fit$lambdaseq), 5)
})

test_that("SPRING stops on negative values", {
  dat <- matrix(c(-1, 1, 2, 3), nrow = 2)
  expect_error(SPRING(dat, nlambda = 3, rep.num = 3, verbose = FALSE))
})

test_that("SPRING warns when quantitative data looks normalized", {
  data(QMP)
  RMP <- QMP[, 1:10] / rowSums(QMP[, 1:10])
  expect_warning(
    SPRING(RMP, quantitative = TRUE, nlambda = 5, rep.num = 5, verbose = FALSE),
    regexp = "normalized"
  )
})

test_that("SPRING is reproducible with the same seed", {
  data(QMP)
  set.seed(42)
  fit1 <- spring_small()
  set.seed(42)
  fit2 <- spring_small()
  expect_equal(fit1$output$stars$opt.index, fit2$output$stars$opt.index)
})
