library(testthat)
library(SCassist)

test_that("SCassist_analyze_quality handles errors", {
  expect_error(SCassist_analyze_quality("nonexistent_object")) 
})