library(testthat)
library(SCassist)

test_that("SCassist_recommend_k handles errors", {
  expect_error(SCassist_recommend_k("nonexistent_object", num_pcs = 6)) 
})