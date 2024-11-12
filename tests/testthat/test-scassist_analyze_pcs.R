library(testthat)
library(SCassist)

test_that("SCassist_analyze_pcs handles errors", {
  expect_error(SCassist_analyze_pcs("nonexistent_object")) 
})