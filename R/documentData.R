#' A test DESeqDataSet
#'
#' test data to check the functions
#'
#' @usage data(testData)
#' @format a  DESeqDataSet
#' @return a  DESeqDataSet
"testData"

#' An example sample table data.frame
#'
#' An example sample table
#'
#' @usage data(exampleSampleTable)
#'
#' @format a data.frame
#' @return a data.frame
#' @seealso   \link{exampleChipData}   \link{exampleInputData}
"exampleSampleTable"

#' An example input data.
#'
#' An example input data set that can be used with \code{\link{importDataFromMatrices}}.
#' Its associated sample table  can be accessed via data(data(exampleSampleTable).
#'
#'
#' @usage data(exampleInputData)
#'
#' @format a matrix
#' @return a matrix
#' @seealso  \link{exampleSampleTable}  \link{exampleChipData}
"exampleInputData"


#' An example ChIP data.
#'
#' An example Chip data set that can be used with \code{\link{importDataFromMatrices}}.
#' Its associated sample table can be accessed via data(data(exampleSampleTable).
#'
#' @usage data(exampleChipData)
#'
#' @format a matrix
#' @return a matrix
#' @seealso  \link{exampleSampleTable}  \link{exampleInputData}
"exampleChipData"

