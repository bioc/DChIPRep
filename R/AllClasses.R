#' @import DESeq2
#' @rdname DChIPRepResults
#' @export

# set the class and the validity check
setClass("DChIPRepResults",
    slots = c(DESeq2Data = "DESeqDataSet",
    FDRresults = "list",
    results = "data.frame"))

setValidity( "DChIPRepResults", function( object ) {
    if (class(DESeq2Data(object)) != "DESeqDataSet" ) {
        return( "The  data is not a DESeqDataSet!" )
        }

    if (!validObject(DESeq2Data(object))  ) {
        return( "The  data is not a valid DESeqDataSet!" )
        }
    })


#' DChIPRepResults object and constructor
#'
#' The \code{DChIPRepResults} contains a DESeqDataSet as obtained after the
#' initial
#' import.
#'
#' @param object A DESeqDataSet
#'
#' @return A DChIPRepResult object.
#'
#'
#'
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#'
#' @rdname DChIPRepResults
#' @import DESeq2
#' @export

### constructor

DChIPRepResults <- function(object) {

    new("DChIPRepResults",
        DESeq2Data = object,
        results = data.frame(),
        FDRresults = list()
)

}

