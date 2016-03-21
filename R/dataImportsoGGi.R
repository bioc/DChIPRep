#' Import the data from bam files directly
#'
#' This function imports the data from .bam files directly. It will
#' return a matrix with one column per .bam file and the respective counts per
#' postion in the rows. It uses the function \code{\link[soGGi]{regionPlot}} 
#' from the package
#' \code{soGGi}.
#'
#' In the example below, we use a subsampled .bam file (0.1 \% of the reads)
#' from the Galonska et. al. WCE (whole cell extract) H3Kme3 data and associated
#' TSS near identified peaks. For additional details on the data, 
#' see \link{input_galonska}
#' and \link{TSS_galonska}.
#'
#' @include AllGenerics.R
#'
#'
#' @param bam_paths a character vector of paths to the bam file(s) 
#' to be imported.
#' @param TSS a GRanges (\link[GenomicRanges]{GenomicRanges-class}) (or a class
#' that inherets from it)
#' object containing the TSS of interest.
#' @param fragment_lengths an integer vector of fragment lengths,
#' @param sample_ids a character vector of sample ids for the .bam files.
#' This can also be a factor.
#' @param distanceUp Distance upstream from centre of the TSS provided.
#' @param distanceDown Distance downstream from centre of the TSS provided.
#'
#' @param ... additional arguments passed to \link[soGGi]{regionPlot}.
#'
#'
#' @return a matrix that contains the postion-wise profiles per .bam 
#' file in the colmuns.
#'
#'
#' @importFrom soGGi regionPlot
#' @importFrom purrr map pmap map_lgl
#'
#'
#' @export
#'
#' @seealso \code{\link[soGGi]{regionPlot}} \link{input_galonska} \link{TSS_galonska} \link{sample_table_galonska}
#'
#'
#' @examples
#' \dontrun{
#' data(sample_table_galonska)
#' data(TSS_galonska)
#' bam_dir <- file.path(system.file("extdata", package="DChIPRep"))
#' wce_bam <- "subsampled_0001_pc_SRR2144628_WCE_bowtie2_mapped-only_XS-filt_no-dups.bam"
#' mat_wce <- importData_soGGi(bam_paths = file.path(bam_dir, wce_bam),
#'                            TSS = TSS_galonska,
#'                            fragment_lengths = sample_table_galonska$input_fragment_length[1],
#'                            sample_ids =  sample_table_galonska$input[1],
#'                            paired = FALSE,
#'                            removeDup=FALSE
#' )
#' head(mat_wce)
#' }
#'



importData_soGGi <- function(bam_paths, TSS, fragment_lengths, sample_ids,
                            distanceUp = 1000, distanceDown = 1500, ...){

    if( !all(map_lgl(bam_paths, is.character)) ){
        stop("The file paths needs to be given as a character vector")
    }

    if( !all(file.exists(bam_paths)) ){
        stop("One or more of the specified .bam files do not exist")
    }

    if( !(inherits(TSS, "GenomicRanges")) ){
        stop("TSS has to inherit from GenomicRanges")
    }

    if( !all(map_lgl(fragment_lengths, is.count)) ){
        stop("fragment_lengths need to be integers")
    }


    if( !( all(map_lgl(sample_ids, is.character)) | all(map_lgl(sample_ids, is.factor)) )  ){
        stop("sample_ids have to be characters")
    }


    .soGGi_import_func <- function(bam, fragment_length, sample_id){

        regionPlot(bamFile = bam, testRanges = TSS, style="point",
                    distanceUp = distanceUp, distanceDown = distanceDown,
                    format="bam", FragmentLength = fragment_length,
                    samplename = sample_id)

    }


## create list of args to iterate over

args_counting <- list(bam_paths, fragment_lengths, sample_ids)
names(args_counting) <- c("bam", "fragment_length", "sample_id")

count_tables <- pmap(args_counting, .soGGi_import_func)

names(count_tables) <- args_counting$sample_id

# function to extract the data matrix and summarize it across positions

.get_data_matrix <- function(count_tables){

    ret <-  map(count_tables, function(X){

        mat <- assay(X)
        summarizeCountsPerPosition(mat)

    })

    ret <- do.call(cbind, ret)

    }

mat <- .get_data_matrix(count_tables)

}


