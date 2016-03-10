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
#' An example ChIP data set that can be used with \code{\link{importDataFromMatrices}}.
#' Its associated sample table can be accessed via data(data(exampleSampleTable).
#'
#' @usage data(exampleChipData)
#'
#' @format a matrix
#' @return a matrix
#' @seealso  \link{exampleSampleTable}  \link{exampleInputData}
"exampleChipData"


#' TSS around called peak regions  from Galonska et. al.
#'
#' This data contains mouse mm9 TSS close to called peak regions for H3Kme3 from Galonksa
#' et al. The original peak lists are from are from GEO (GSE56312) and have been
#' merged into a common peaklist and then annotated to the closest mm9 TSS using
#' \link[ChIPpeakAnno]{annotatePeakInBatch}.
#'
#' @usage data(exampleChipData)
#'
#' @import ChIPpeakAnno
#'
#' @format an annoGR object from the package \code{ChIPpeakAnno}.
#' @return an annoGR object from the package \code{ChIPpeakAnno}.
#' @seealso  \link{chip_galonska} \link{input_galonska}  \link{sample_table_galonska}
"TSS_galonska"



#' Another example ChIP data set that can be used with \code{\link{importDataFromMatrices}} from rom Galonska et. al., 2015.
#'
#' The containes the H3Kme3 data in the serum and 24h_2i conditions from Galonska et. al., 2015
#' as well as the whole cell extract data, which is treated as input for all
#' four samples.
#' The data were downloaded  from the SRA at the european
#' nucleotide archive (ENA, accession PRJNA242892).
#' The reads were aligned ot the mm9 reference genome using bowtie2
#' (Langmead and Salzberg, 2012) with
#' default options. Then, filtering of unmapped, low mapping quality (< 10), duplicated and
#' multi-mapping reads was performed with Picard tools.
#' The fragment length was inferred using
#' cross correlation plots from SPP (Kharchenko, et. al., 2008).
#'
#' @references
#'
#' Galonska, Christina, Michael J. Ziller, Rahul Karnik, and Alexander Meissner. 2015. "Ground State Conditions Induce Rapid Reorganization of Core Pluripotency Factor Binding Before Global Epigenetic Reprogramming." Cell Stem Cell 17 (4). Elsevier BV: 462-70.
#' \url{http://dx.doi.org/10.1016/j.stem.2015.07.005}.

#'Kharchenko, Peter V, Michael Y Tolstorukov, and Peter J Park. 2008. "Design and Analysis of ChIP-Seq Experiments for DNA-Binding Proteins." Nat Biotechnol 26 (12). Nature Publishing Group: 1351-9.  \url{http://dx.doi.org/10.1038/nbt.1508}.

#'Langmead, Ben, and Steven L Salzberg. 2012. "Fast Gapped-Read Alignment with Bowtie 2." Nature Methods 9 (4). Nature Publishing Group: 357-59. \url{http://dx.doi.org/10.1038/nmeth.1923}.

#'Picard Tools - by Broad Institute. 2016. \url{http://broadinstittue.github.io/picard/}.

#'
#' @usage data(chip_galonska)
#'
#' @format a matrix
#' @return a matrix
#' @seealso  \link{sample_table_galonska}  \link{input_galonska} \link{TSS_galonska}
"chip_galonska"

#' Another example Input data set that can be used with \code{\link{importDataFromMatrices}} from rom Galonska et. al., 2015.
#'
#' The matrix containes the whole cell extract (WCE) data for H3Kme3 from the paper
#' in each of the four columns, since
#' this is the only inpu data provided for all 4 samples. For additional
#'  information see the documentation
#'  of \link{chip_galonska}.
#'
#' @usage data(input_galonska)
#'
#' @format a matrix
#' @return a matrix
#' @seealso  \link{chip_galonska}  \link{sample_table_galonska} \link{TSS_galonska}
"input_galonska"

#' Another example sample table based on data  from rom Galonska et. al., 2015.
#'
#' This table contains the sample annotation for the H3Kme3 data from
#'  Galonska et. al., 2015. For additional information see the documentation
#'  of \link{chip_galonska}.
#'
#' @usage data(sample_table_galonska)
#'
#' @format a data.frame
#' @return a data.frame
#' @seealso  \link{chip_galonska}  \link{input_galonska} \link{TSS_galonska}
"sample_table_galonska"

