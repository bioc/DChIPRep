library(DESeq2)
#library(ggplot2)

context("DChIPRep plotting functions")

test_that("plottingFunctions_work", {
    data(testData)
    dcr <- DChIPRepResults(testData)
    dcr <- runTesting(dcr)

    resS <-plotSignificance(dcr)
    resP <-plotProfiles(dcr)

    expect_is(resS, "ggplot")
    expect_is(resP, "ggplot")
})
