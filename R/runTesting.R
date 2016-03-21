runTesting.DChIPRepResults <- function(object,
                                        lfcThreshold = 0.05,
                                        plotFDR = FALSE,...){
    # Dispersion estimation
    DESeq2Data <- estimateDispersions(DESeq2Data(object), fitType = "local")

    # Stat. test
    DESeq2Data <- nbinomWaldTest(DESeq2Data)

    ### Get the results
    res <- as.data.frame(DESeq2::results(DESeq2Data,lfcThreshold = lfcThreshold,
                            altHypothesis="greaterAbs",
                            independentFiltering = FALSE,
                            format = "DataFrame", ...))



# FDR computations

    FDR <- fdrtool(res$pvalue, statistic = "pvalue", plot = plotFDR)

    res$lfdr <- FDR$lfdr


    # add FDR to the results and to the object

    FDRresults(object) <- FDR
    resultsDChIPRep(object) <- res
    return(object)
}


#' Run the tests on a DChIPRepResults object.
#'
#' This function runs the testing on a \code{DChIPRepResults} object. It adds
#' the FDR calculations and the result table to the \code{DChIPRepResults}
#' object.
#'
#' @include AllGenerics.R
#' @name runTesting
#' @rdname runTesting
#' @aliases runTesting runTesting,DChIPRepResults-method
#' @param object A \link{DChIPRepResults} object.
#' @param lfcThreshold A non-negative threshold value, which determines the null
#' hypothesis. The null hypothesis is
#' \eqn{H_0: |log2(FC)| > } \code{lfcThreshold}
#'
#' @param plotFDR If set to TRUE a plot showing the
#' estimated FDRs will be displayed
#' @param ... not used currently
#'
#'
#' @return a modified \link{DChIPRepResults} object containing the 
#' testing results
#' @seealso \code{\link{resultsDChIPRep}}
#' @export
#' @import DESeq2
#' @importFrom fdrtool fdrtool
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' dcr <- runTesting(dcr)
setMethod("runTesting", signature(object="DChIPRepResults"),
            runTesting.DChIPRepResults)



