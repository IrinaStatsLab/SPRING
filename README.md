
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

The faster version of latent correlation computation part is now fully
available and implemented to the R package `SPRING`. 
## Installation

``` r
# install.packages("devtools")
devtools::install_github("GraceYoon/SPRING")
```

## Example

``` r
library(SPRING)
data("QMP") # load the data available from this package, containing 106 samples and 91 OTUs.

# Apply SPRING on QMP data.
fit.spring <- SPRING(QMP, Rmethod = "approx", quantitative = TRUE, 
                     lambdaseq = "data-specific", nlambda = 50, rep.num = 50)
# With Rmethod = "original", this takes around 23 minutes.
# With Rmethod = "approx", this takes around 2.23 minutes. 
# More details on the comparison of accuracy and speed ("original" vs. "approx")
# are available on the above arXiv reference.

# StARS-selected lambda index based on the threshold (default = 0.01)
opt.K <- fit.spring$output$stars$opt.index
# Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
adj.K <- as.matrix(fit.spring$fit$est$path[[opt.K]])
# Estimated partial correlation coefficient, same as negative precision matrix.
pcor.K <- as.matrix(SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K]], mode = 'maxabs'))
```
