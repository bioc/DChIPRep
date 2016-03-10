
context("soGGi based data import")

# comment this since it takes too long ...
# test_that("soGGi_based_import works", {
#   data(sample_table_galonska)
#   data(TSS_galonska)
#   bam_dir <- file.path(system.file("extdata", package="DChIPRep"))
#   wce_bam <- "subsampled_0001_pc_SRR2144628_WCE_bowtie2_mapped-only_XS-filt_no-dups.bam"
#   mat_wce <- importData_soGGi(bam_paths = file.path(bam_dir, wce_bam),
#                               TSS = TSS_galonska,
#                               fragment_lengths = sample_table_galonska$input_fragment_length[1],
#                               sample_ids =  sample_table_galonska$input[1],
#                            paired = FALSE,
#                              removeDup=FALSE
#   )
#  expect_is(mat_wce, "matrix")
# })


test_that("bam_paths checks work", {
    data(sample_table_galonska)
    data(TSS_galonska)
    bam_dir <- file.path(system.file("extdata", package="DChIPRep"))
    wce_bam <- "subsampled_0001_pc_SRR2144628_WCE_bowtie2_mapped-only_XS-filt_no-dups.bam"

    expect_error(mat_wce <- importData_soGGi(bam_paths = 0,
                                TSS = TSS_galonska,
                                fragment_lengths = sample_table_galonska$input_fragment_length[1],
                                sample_ids =  sample_table_galonska$input[1],
                                paired = FALSE,
                                removeDup=FALSE),
                 "The file paths needs to be given as a character vector")

    expect_error(mat_wce <- importData_soGGi(bam_paths = "path/does/not/exist",
                                             TSS = TSS_galonska,
                                             fragment_lengths = sample_table_galonska$input_fragment_length[1],
                                             sample_ids =  sample_table_galonska$input[1],
                                             paired = FALSE,
                                             removeDup=FALSE),
                 "One or more of the specified .bam files do not exist")
})


test_that("other input data checks work", {
    data(sample_table_galonska)
    data(TSS_galonska)
    bam_dir <- file.path(system.file("extdata", package="DChIPRep"))
    wce_bam <- "subsampled_0001_pc_SRR2144628_WCE_bowtie2_mapped-only_XS-filt_no-dups.bam"

    expect_error(mat_wce <- importData_soGGi(bam_paths = file.path(bam_dir, wce_bam),
                                             TSS = c(1,2,3),
                                             fragment_lengths = sample_table_galonska$input_fragment_length[1],
                                             sample_ids =  sample_table_galonska$input[1],
                                             paired = FALSE,
                                             removeDup=FALSE),
                 "TSS has to inherit from GenomicRanges")

    expect_error(mat_wce <- importData_soGGi(bam_paths = file.path(bam_dir, wce_bam),
                                             TSS = TSS_galonska,
                                             fragment_lengths = c(-10, 150),
                                             sample_ids =  sample_table_galonska$input[1],
                                             paired = FALSE,
                                             removeDup=FALSE),
                 "fragment_lengths need to be integers")

    expect_error(mat_wce <- importData_soGGi(bam_paths = file.path(bam_dir, wce_bam),
                                             TSS = TSS_galonska,
                                             fragment_lengths = sample_table_galonska$input_fragment_length[1],
                                             sample_ids =  c(55, 77),
                                             paired = FALSE,
                                             removeDup=FALSE),
                 "sample_ids have to be characters")
})