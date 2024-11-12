library(testthat)
library(SCassist)

test_that("SCassist_analyze_and_annotate handles errors for non-existent object", {
  expect_error(SCassist_analyze_and_annotate(mock_all_markers, seurat_object_name = "nonexistent_object"))
}) 