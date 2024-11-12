library(testthat)
library(SCassist)

test_that("SCassist_recommend_res handles errors", {
  expect_error(SCassist_recommend_res("nonexistent_object")) 
})