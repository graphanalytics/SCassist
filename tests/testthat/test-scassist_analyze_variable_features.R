library(testthat)
library(SCassist)

test_that("SCassist_analyze_variable_features handles errors", {
  expect_error(SCassist_analyze_variable_features("nonexistent_object")) 
})