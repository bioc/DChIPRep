library(DESeq2)

context("DChIPRep data import")

test_that("dataImport_works", {
    directory <- file.path(system.file("extdata", package="DChIPRep"))
    data(exampleSampleTable)
    res <- importData(exampleSampleTable, directory)
    expect_is(res, "DChIPRepResults")
})

