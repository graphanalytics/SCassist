library(testthat)
library(SCassist)

test_that("SCassist_recommend_normalization handles errors", {
  expect_error(SCassist_recommend_normalization("nonexistent_object")) 
})