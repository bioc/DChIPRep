context("data import tests")

test_that("dataImportColumnTests_work", {
    data(exampleSampleTable)
    data(exampleInputData)
    data(exampleSampleTable)
    directory <- file.path(system.file("extdata", package="DChIPRep"))
    exampleSampleTable$ChIP <- NULL
    expect_error(importedData <- importData(exampleSampleTable, directory),
    "The sample table must contain the columns Input and ChIP !")
})


test_that("dataImportFileExistsTests_work", {
    data(exampleSampleTable)
    data(exampleInputData)
    data(exampleSampleTable)
    directory <- file.path(system.file("extdata", package="DChIPRep"))
    exampleSampleTable$Input <- rep("XXX", dim(exampleSampleTable)[1] )
    expect_error(importedData <- importData(exampleSampleTable, directory),
    "Some or all Input/ChIP files do not exist!")
})



test_that("importDataFromMatricesInputFormatCheck_works", {
    data(exampleSampleTable)
    data(exampleInputData)
    data(exampleChipData)

    exampleSampleTable <- numeric(0)

    expect_error(imDataFromMatrices <- importDataFromMatrices(
    inputData = exampleInputData,
    chipData = exampleChipData,
    sampleTable = exampleSampleTable),
    "Input sample table is not a data frame")
})

test_that("importDataFromMatricesColumnCheck_works", {
    data(exampleSampleTable)
    data(exampleInputData)
    data(exampleChipData)

    exampleSampleTable$upstream <- NULL
    exampleSampleTable$sampleID <- NULL

    expect_error(imDataFromMatrices <- importDataFromMatrices(
    inputData = exampleInputData,
    chipData = exampleChipData,
    sampleTable = exampleSampleTable),
    "The sample table needs to contain columns named, \n upstream, downstream, condition and sampleID")
})

test_that("importDataFromMatricesUpstreamDownstreamCheck_works", {
    data(exampleSampleTable)
    data(exampleInputData)
    data(exampleChipData)
    exampleSampleTable$upstream <- rep(-4444, dim(exampleSampleTable)[1])
    exampleSampleTable$downstream <- 77.7


    expect_error(imDataFromMatrices <- importDataFromMatrices(
    inputData = exampleInputData,
    chipData = exampleChipData,
    sampleTable = exampleSampleTable),
    "Upstream and/or downstream positions are not unique \n and/or  are not positive integers")
})

