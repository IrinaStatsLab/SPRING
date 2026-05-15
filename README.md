
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPRING

<!-- badges: start -->

<!-- badges: end -->

The R package `SPRING` (Semi-Parametric Rank-based approach for
INference in Graphical model) estimates sparse microbial association
networks using rank-based correlation with sparse graphical modeling
techniques. The corresponding reference is

Yoon G., Gaynanova I. and Müller C.L. (2019) [Microbial Networks in
SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation
Estimation for Quantitative Microbiome
Data](https://www.frontiersin.org/articles/10.3389/fgene.2019.00516/full).
*Frontiers in Genetics*, 10:516.

Latent correlation estimation is provided by the R package `latentcor`
(Huang Z., Guo M., Müller C.L. and Gaynanova I. (2021) [latentcor: An R
Package for Estimating Latent Correlations from Mixed Data
Types](https://doi.org/10.21105/joss.03634). *Journal of Open Source
Software*, 6(65), 3634).

## Installation

The package is being submitted to Bioconductor. Until it is available
there, install from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("IrinaStatsLab/SPRING")
```

## Example

``` r
library(SPRING)
data("QMP") # 106 samples and 91 taxa

# Apply SPRING on QMP data
fit <- SPRING(QMP, quantitative = TRUE, nlambda = 50, rep.num = 50,
              verbose = FALSE)

# StARS-selected optimal lambda index (default threshold = 0.1)
opt.index <- fit$output$stars$opt.index

# Estimated adjacency matrix (1 = edge, 0 = no edge)
adj <- as.matrix(fit$fit$refit$stars)

# Symmetrized MB coefficients at optimal lambda
beta <- as.matrix(fit$fit$est$beta[[opt.index]])
beta_sym <- (beta + t(beta)) / 2
```
