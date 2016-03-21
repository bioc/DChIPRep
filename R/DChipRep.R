#' @import methods GenomicRanges ggplot2
NULL

#' DChIPRep: A package for differential analysis of histone modification ChIP-Seq profiles
#'
#' The DChIPRep package provides functions to perform a differential analysis
#' of histone modification profiles at base-pair resolution
#'
#' @section DChIPRep functions:
#' The \code{DChIPRep} packages provides functions for data import
#' \code{\link{importData}}
#' and performing
#' position wise tests. After data import, a \code{\link{DChIPRepResults}}
#' object on which the function \code{\link{runTesting}}
#' is run to perform the tests and add the result to the object. Then, plots
#' can be created from this object. See the vignette for additional details:
#' \code{vignette("DChIPRepVignette")}
#'
#'
#' @docType package
#' @name DChIPRep
NULL
#> NULL
