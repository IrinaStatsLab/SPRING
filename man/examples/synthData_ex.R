data(QMP)

# Fast example: small subset (5 taxa), identity covariance
Sigma_small <- diag(5)
X_small <- synthData_from_ecdf(QMP[, 1:5], Sigma = Sigma_small, n = 20,
                               seed = 10010)

# Full example with a prescribed graph structure -- slow on all 91 taxa
\donttest{
set.seed(12345)
p <- ncol(QMP)
graph <- SpiecEasi::make_graph("cluster", p, 2 * p)
Prec  <- SpiecEasi::graph2prec(graph)
Sigma <- cov2cor(SpiecEasi::prec2cov(Prec))
X <- synthData_from_ecdf(QMP, Sigma = Sigma, n = 100, seed = 10010)
}
