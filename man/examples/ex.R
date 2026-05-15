data(QMP)

# Fast example: small subset, few subsamples
fit <- SPRING(QMP[, 1:10], quantitative = TRUE, nlambda = 10, rep.num = 10,
              verbose = FALSE)
