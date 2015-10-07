
context("DChIPRep testing test")

test_that("runTesting_works", {
  data(testData)
  dcr <- DChIPRepResults(testData)
  dcr <- runTesting(dcr)
  expect_is(resultsDChIPRep(dcr), "data.frame")
  expect_is(FDRresults(dcr),"list")
})

