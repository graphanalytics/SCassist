library(testthat)
library(SCassist)

test_that("SCassist_analyze_enrichment handles errors for unsupported organism", {
  expect_error(SCassist_analyze_enrichment(organism = "invalid_organism", markers = mock_markers))
})