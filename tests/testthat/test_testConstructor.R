library(S4Vectors)
library(DESeq2)

context("DChIPRepResults constructor and accesors testing")

test_that("constructor_works", {
  countData <- matrix(1:100,ncol=4)
  condition <- factor(c("A","A","B","B"))
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
  res <- DChIPRepResults(dds)
  expect_true(class(res) ==  "DChIPRepResults" )
})


test_that("data_accessor_works", {
  countData <- matrix(1:100,ncol=4)
  condition <- factor(c("A","A","B","B"))
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
  res <- DChIPRepResults(dds)
  expect_true(class(DESeq2Data(res)) ==  "DESeqDataSet" )
})


test_that("show_method_works", {
  countData <- matrix(1:100,ncol=4)
  condition <- factor(c("A","A","B","B"))
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
  dds <- estimateSizeFactors(dds)
  res <- DChIPRepResults(dds)
  res <- runTesting(res)
  res
})
