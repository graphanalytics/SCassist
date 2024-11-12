library(testthat)
library(SCassist)

test_that("SCassist_recommend_pcs handles errors", {
  expect_error(SCassist_recommend_pcs("nonexistent_object")) 
})