data(QMP)

# Fast example: small subset, few subsamples
fit <- SPRING(QMP[, 1:10], quantitative = TRUE, nlambda = 10, rep.num = 10,
              verbose = FALSE)

# Compositional data: apply mclr transformation internally
\donttest{
RMP <- QMP / rowSums(QMP)
fit_comp <- SPRING(RMP, quantitative = FALSE, nlambda = 10, rep.num = 10,
                   verbose = FALSE)
}

# Full analysis as in Yoon et al. (2019) -- slow, ~20 min on full data
\donttest{
fit_full <- SPRING(QMP, quantitative = TRUE, lambdaseq = "data-specific",
                   nlambda = 50, rep.num = 50, verbose = FALSE)
}
