data(QMP)

# Fast example: small subset (5 taxa), identity covariance
set.seed(10010)
Sigma_small <- diag(5)
X_small <- synthData_from_ecdf(QMP[, 1:5], Sigma = Sigma_small, n = 20)
