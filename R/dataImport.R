#' Helper function to turn a  data.frame  into a matrix
#' and remove the ID column.
#'
#' This function takes a data.frame, with the genomic features
#' (e.g. transcripts or genes) in the rows and the positions upstream and
#' downstream of
#' the TSS in the columns as well as a column ID containing a genomic feature ID
#' and returns the data.frame with the ID column removed. The input for this
#' function
#' are tables obtained after running the Python import script.
#'
#'
#' @param  df the input data frame with positions in the columns and the genomic
#' features in the rows.
#' @param ID  the name of the ID column to be removed.
#'
#'
#'
#' @export
#' @return a matrix with the ID column removed
#' @examples
#' data(exampleSampleTable)
#' directory <- file.path(system.file("extdata", package="DChIPRep"))
#' df  <- lapply(file.path(directory, exampleSampleTable$Input), 
#' read.delim)[[1]]
#' mat <- getMATfromDataFrame(df)


getMATfromDataFrame <- function(df, ID="name"){

  if(!(ID %in% names(df))){
  stop(paste0("There is no column called ", ID, " in \n the input data"))
  }

 df[,ID] <- NULL
 as.matrix(df)
}




#' Helper function to summarize the counts per position
#'
#' This function takes a matrix of counts, with the genomic features
#' (e.g. transcripts or genes) in the rows  and the positions upstream
#' and downstream of
#' the TSS in the columns and returns a vector with the summarized counts per
#' position.
#'
#' The summary per condition is computed as a trimmed mean per position. First,
#' counts greater than \code{ct}  are removed and then a trimmed mean with
#' a trimming
#' percentage of \code{trim} is computed on the log scale. This mean is
#' then exponentiated
#' again and multiplied by the total number of features passing the threshold
#' \code{ct}
#' per position. If a position contains only zero counts, its mean
#' is returned as zero.
#'
#' @importFrom plyr llply
#'
#' @param  mat the input matrix with positions in the columns and the genomic
#' features in the rows.
#' @param ct the count threshold to use.
#' @param trim  the trimming percentage for the trimmed mean.
#'
#' @export
#' @return a vector containing the summarized counts per condition
#'
#'
#' @examples
#' data(exampleSampleTable)
#' directory <- file.path(system.file("extdata", package="DChIPRep"))
#' df  <- lapply(file.path(directory, exampleSampleTable$Input), 
#' read.delim)[[1]]
#' mat <- getMATfromDataFrame(df)
#' summaryPerPos <- summarizeCountsPerPosition(mat)


summarizeCountsPerPosition <- function(mat, ct=0, trim=0.15){
  # Get log counts
  tp <- log(mat + 1)

  # Remove 0 counts and counts below ct
  tp[,colMeans(abs(tp)) < .Machine$double.eps ^ 0.5] <- NA
  tp[tp <= ct] <- NA

  # Get the number of genes with counts > 0 at each position
  factors <- colSums(!is.na(tp))

  # Get the mean counts at each position
  matMEANS <- vapply(split(tp, col(tp)), mean, double(1L), trim = trim, 
  na.rm = TRUE)
  matMEANS[is.na(matMEANS)] <- 0

  # Multiply by the number of genes with counts > 0
  ceiling(factors * expm1(matMEANS))
}





# Function to run tests on the sample table
testSampleTable <- function(sampleTable){

   # test whether the sample Table is a data frame
  if( !"data.frame" %in% class(sampleTable)){
    stop("Input sample table is not a data frame")
  }

   # test whether the sample table contains neccessary columns
  testCI <- function(sampleTable) {
    all("sampleID" %in% names(sampleTable),
    "condition" %in% names(sampleTable),
    "upstream" %in% names(sampleTable),
    "downstream" %in% names(sampleTable))
  }
  on_failure(testCI) <- function(call, env) {
    paste0("The sample table needs to contain columns named, \n upstream, downstream, condition and sampleID")
  }
   assert_that(testCI(sampleTable))

    # test uniqueness and equality of the upstream and downstream vectors
  testUniqueness <- function(up, down) {
    all(sapply(up, is.count),
    sapply(down, is.count),
    length(unique(up)) == 1,
    length(unique(down)) == 1)
  }
  on_failure(testUniqueness) <- function(call, env) {
    paste0("Upstream and/or downstream positions are not unique \n and/or  are not positive integers")
  }
  assert_that(testUniqueness(sampleTable$upstream,
                sampleTable$downstream))

}


#' Import the data from ChiP and input matrices
#'
#' This function imports the data from two matrices that
#' contain counts summarized
#' per position. It computes the normalization factors from the input
#' (one per position)  and creates a DChIPRepResults
#' object.
#'
#'
#' The normalization factors are computed as
#' \code{t(t(inputData) * (covC/covI)) },  Where covC and covI contain the total
#' sum of the ChIP and the input samples. Zero normalization factors can arise
#' if the input has zero counts for certain
#' positions. That's why input values equal to zero are set to 1
#' in order to always
#' obtain valid normalizationFactors.
#'
#' @include AllGenerics.R
#'
#'
#' @param inputData a matrix containing the counts for the input per position.
#' @param chipData a matrix containing the counts for the ChIP per position.
#' @param sampleTable a data.frame that has to contain the columns sampleID,
#' upstream, downstream and condition. Each row of the table describes one
#' experimental sample. See \code{data(exampleSampleTable)} for an example 
#' table.
#' and the vignette for further information.
#'
#'
#'
#' @return a DChIPRepResults object containing the imported data as a
#'\code{\link[DESeq2]{DESeqDataSet}}.
#'
#'
#' @import DESeq2
#'
#' @export
#' @examples
#' data(exampleSampleTable)
#' data(exampleInputData)
#' data(exampleChipData)
#' imDataFromMatrices <- importDataFromMatrices(inputData = exampleInputData,
#' chipData = exampleChipData,
#' sampleTable = exampleSampleTable)

importDataFromMatrices <- function(inputData, chipData, sampleTable){
    testSampleTable(sampleTable)

    expGroups <- sampleTable$condition

    # input and ChIP coverage
    covI <- apply(inputData, 2, function(x){ tmp = as.numeric(x);  sum(tmp)})
    covC <- apply(chipData, 2, function(x){ tmp = as.numeric(x);  sum(tmp)})

    # replace zeros with 1s in the input  data

    inputData[ abs(inputData) < .Machine$double.eps ^ 0.5 ] <- 1

    # multiply input ROWS with ChIP/Input ratio vector


    normFactors <- t(t(inputData) * (covC/covI))


    gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    # normalize to geometric mean 1
    normFactors <- normFactors/gm_mean(normFactors)

    mode(chipData) <- "integer"

    # turn grouping vector into a factor
    sampleTable$condition  <- as.factor(sampleTable$condition)
    condition <- as.factor(sampleTable$condition)


    if(!all.equal(length(levels(condition)),2) ){
    stop("There must be exactly 2 experimental groups !")
    }
    
    # check whether the colnames of the input matrices match the sample IDs
    if(any(sampleTable$sampleID != colnames(inputData)) &  any(sampleTable$sampleID != colnames(chipData)) ){
      stop("The column names of your ChIP and Input matrices have to be equal to the sampleID column")
    }

    ## get the positions as row names
    rownames(inputData) <- paste0("Pos_",seq(from =-unique(sampleTable$upstream)
                            , to = unique(sampleTable$downstream), by = 1))
    rownames(chipData) <- paste0("Pos_", seq(from =-unique(sampleTable$upstream)
                            , to = unique(sampleTable$downstream), by = 1))


    ## check whether the sample IDs contain duplicates
    if(all(duplicated(sampleTable$sampleID))){

    stop("The sample IDs are not unique.")

    }
    ## get the sample names as column names
    colnames(inputData) <- sampleTable$sampleID
    colnames(chipData) <- sampleTable$sampleID


    # create a DESeq2 object
    DESeq2Data <- DESeqDataSetFromMatrix(countData = chipData,
                                    colData = data.frame(sampleTable),
                                    design = formula(~ condition))
    normalizationFactors(DESeq2Data) <- normFactors

    DChIPRepResults(DESeq2Data)
}

#' Import the data after running the Python script
#'
#' This function imports the data from the count table files as returned by
#' the accompanying Python script.
#'
#'
#' @include AllGenerics.R
#'
#' @importFrom plyr ldply
#' @importFrom stats aggregate.data.frame formula
#' @importFrom utils read.delim
#'
#' @import assertthat
#'
#' @param sampleTable a data.frame that has to contain the columns ChiP, Input,
#' sampleID,
#' upstream, downstream and condition. Each row of the table describes one
#' experimental sample. Each row of the table describes one
#' experimental sample. See \code{data(exampleSampleTable)} for an example 
#' table.
#' and the vignette for further information.
#'
#' @param directory  the directory relative to which the filenames are specified
#' given as a character.
#' @param ID character giving the name of the feature identifier
#' column in the count tables.
#'
#' Defaults to \code{"name"}
#' @param ... parameters passed to \code{\link{summarizeCountsPerPosition}}
#'
#' @return a \code{DChIPRepResults} object containg the imported data as a
#' \code{\link[DESeq2]{DESeqDataSet}}.
#'
#' @export
#' @examples
#' data(exampleSampleTable)
#' directory <- file.path(system.file("extdata", package="DChIPRep"))
#' importedData <- importData(exampleSampleTable, directory)

importData <- function(sampleTable, directory="", ID="name", ...){


if(!all("Input" %in% names(sampleTable),
    "ChIP" %in% names(sampleTable))){
    stop("The sample table must contain the columns Input and ChIP !")
    }


stopifnot(is.character(ID))
stopifnot(is.character(directory))


if( !all(file.exists(file.path(directory, sampleTable$Input)),
        file.exists(file.path(directory, sampleTable$ChIP)))){
        stop("Some or all Input/ChIP files do not exist!")
        }


inputData <- lapply(file.path(directory, sampleTable$Input), read.delim)
chipData <- lapply(file.path(directory, sampleTable$ChIP), read.delim)



inputData <- lapply(inputData, getMATfromDataFrame)
chipData <- lapply(chipData, getMATfromDataFrame)


inputData <- t(as.matrix(ldply(inputData, summarizeCountsPerPosition, ...)))
chipData <- t(as.matrix(ldply(chipData, summarizeCountsPerPosition, ...)))




importDataFromMatrices(inputData, chipData, sampleTable)

}



