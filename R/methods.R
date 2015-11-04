# This file contains several (generic) accession methods


# Accessors for the 'DESeq2Data' slot of a DChIPRepResults object.
DESeq2Data.DChIPRepResults <- function(object){
  object@DESeq2Data
}


#' Accessors for the 'DESeq2Data' slot of a \code{DChIPRepResults} object.
#'
#' The slot contains the DESeqDataSet as it is obtained after
#' the initial data import. The DESeqDataSet contains the
#' counts per position and the normalization factors as computed  using
#' the input counts.
#'
#' @name DESeq2Data
#' @rdname DESeq2Data
#' @param object a \code{DChIPRepResults} object
#'
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' DESeq2Data(dcr)
#'
#' @rdname DESeq2Data
#' @importFrom DESeq2 DESeqDataSet
#' @aliases DESeq2Data DESeq2Data<- DESeq2Data,DChIPRepResults-method
#'
#' @return the DESeq2Data object contained in the DChIPRepResults object
#' @export
setMethod("DESeq2Data", signature(object="DChIPRepResults"),
            DESeq2Data.DChIPRepResults)

#' @param value A DESeqDataSet object
#' @rdname DESeq2Data
#' @importFrom DESeq2 DESeqDataSet
#' @export
setMethod("DESeq2Data<-", signature(object="DChIPRepResults",
                                    value="DESeqDataSet"),
            function( object, value ) {
                object@DESeq2Data <- value
                validObject(object)
                object
            })


## create the  show method

#' prints the DESeq2Data slot of the DChIPRepResults object
#'
#' prints the data
#'
#' @rdname show
#'
#' @importFrom S4Vectors DataFrame
#' @name show
#' @param object A DChIPRepResults object
#' @aliases show show,DChIPRepResults-method
#'
#'
#' @export
#' @return A compact representation of the DChIPRepResults object
#'
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' dcr
#' dcr <- runTesting(dcr)
#' dcr

setMethod("show", signature(object="DChIPRepResults"), function(object) {
                cat("Summary of the DESeqDataSet:\n")
                cat("============================\n")
                show(DESeq2Data(object))
                cat("\n\n\n")
                cat("Results:\n")
                cat("============================\n")
                if(isTRUE(all.equal(dim(resultsDChIPRep(object)), c(0,0)))){
                cat("No results available yet, call runTesting first\n\n")
                }else{
                print(DataFrame(resultsDChIPRep(object)))
                }
            })


#' Accessor and setter for the 'FDRresults' slot of a \code{DChIPRepResults}
#' object.
#'
#' The slot contains the results of the FDR estimation as performed within the
#' function \code{\link{runTesting}}. It is the complete output of the
#' \code{\link[fdrtool]{fdrtool}} function.
#'
#' @name FDRresults
#' @rdname FDRresults
#' @param object a \code{DChIPRepResults} object
#'
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' dcr <- runTesting(dcr)
#' str(FDRresults(dcr))
#'
#' @return a list containing the estimated false discovery rates
#'
#' @rdname FDRresults
#' @importFrom DESeq2 DESeqDataSet
#' @aliases FDRresults FDRresults<- FDRresults,DChIPRepResults-method
#' @export
setMethod("FDRresults", signature(object="DChIPRepResults"),
            function(object){
                object@FDRresults
            })

#' @param value A DESeqDataSet object
#' @rdname FDRresults
#' @importFrom DESeq2 DESeqDataSet
#' @export
setMethod("FDRresults<-", signature(object="DChIPRepResults", value="list"),
            function(object, value){
                object@FDRresults <- value
                validObject(object)
                object
            })

#' Accessors and setter for the 'results' slot of a \code{DChIPRepResults} object.
#'
#' The slot contains the results of the position wise tests in a
#' data.frame after runing
#' the function \code{\link{runTesting}}. It is a modified output of the
#' \code{\link[DESeq2]{results}} function of the \code{DESeq2} package.
#'
#'
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' dcr <- runTesting(dcr)
#' head(resultsDChIPRep(dcr))
#'
#' @name resultsDChIPRep
#' @param object a \code{DChIPRepResults} object
#' @rdname resultsDChIPRep
#' @aliases resultsDChIPRep resultsDChIPRep<- resultsDChIPRep,DChIPRepResults-method
#' @export
#' @return a data.frame containing the results of the position wise tests
#'
setMethod("resultsDChIPRep", signature(object="DChIPRepResults"),
        function(object){
            object@results
        })

#' @param value A DESeqDataSet object
#' @rdname resultsDChIPRep
#' @export
setMethod("resultsDChIPRep<-", signature(object="DChIPRepResults",
                                         value="list"),
        function(object, value){
            object@results <- value
            validObject(object)
            object
        })

