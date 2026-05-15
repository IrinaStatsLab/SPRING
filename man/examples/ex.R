# Load the synthetic count data
data("QMP") # n = 1000 and p = 100 synthetic dataset

# SPRING on Synthetic Data, when assuming the data as quantitative counts.
# The same setting used in Yoon et al. (2019) Frontiers in Genetics.
\dontrun{
# This takes around 23 minutes.
fit.spring <- SPRING(QMP, quantitative = TRUE, lambdaseq = "data-specific",
                     nlambda = 50, seed = 10010, ncores = 2, rep.num = 50)
}

# SPRING on Compositional data. Row sums are scaled to 1. Then, mclr-transformation will be applied.
\dontrun{
compoData <- QMP/rowSums(QMP)
fit.spring <- SPRING(compoData, quantitative = FALSE, lambdaseq = "data-specific",
                     nlambda = 10, rep.num = 10)
}
