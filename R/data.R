#' Synthetic count data
#' @name SynthData
#' @description  SynthData and SynthData2 were generated using empirical cdf of American Gut Project Data. SynthData has scale_free-type-graph structure of size 500 rows and 200 columns, and SynthData2 has cluster-type-graph structure of size 1000 rows and 100 columns.
#'
#' @source
#'
#' Yoon, Gaynanova and Müller (2019) Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data. \emph{Frontiers in Genetics.} 10:516. \url{doi:10.3389/fgene.2019.00516}
#' @format \code{SynthData} is an object of class \code{matrix} with 500 rows and 200 columns. \code{SynthData2} is an object of class \code{matrix} with 1000 rows and 100 columns.
"SynthData"



#' @name SynthData
#' @aliases SynthData2
"SynthData2"

#' Quantitative Microbiome Project data
#'
#' @description Quantitative microbiome count data from Vandeputte et al. (2017)
#'   with 106 healthy subjects (rows) and 91 taxa (columns). Taxa present in
#'   fewer than 30\% of samples are excluded, and only healthy subjects from the
#'   Study cohort and Disease cohort are retained. Column names are taxonomy
#'   labels: genus name where classified, otherwise the most specific named rank
#'   (family, order, class, or phylum), with numeric suffixes disambiguating
#'   multiple unclassified taxa within the same rank (e.g.,
#'   \code{Ruminococcaceae} and \code{Ruminococcaceae_1}).
#'
#' @source
#'
#' Yoon, Gaynanova and Mueller (2019) Microbial Networks in SPRING -
#' Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for
#' Quantitative Microbiome Data. \emph{Frontiers in Genetics}. 10:516.
#' \url{doi:10.3389/fgene.2019.00516}
#'
#' Vandeputte et al. (2017) Quantitative microbiome profiling links gut
#' community variation to microbial load. \emph{Nature}. 551: 507-511.
#' \url{doi:10.1038/nature24460}
#'
"QMP"
