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


test_that("plottingFunctions_work_without_replicates", {
  data(testData)
  dcr <- DChIPRepResults(testData[, 3:4])
  
  # since version 1.22 one has to change the design explicitly
  dds <- DESeq2Data(dcr)
  design(dds) <- ~ 1
  DESeq2Data(dcr) <- dds
  
  dcr <- runTesting(dcr)
  
  resS <-plotSignificance(dcr)
  resP <-plotProfiles(dcr)
  
  expect_is(resS, "ggplot")
  expect_is(resP, "ggplot")
})