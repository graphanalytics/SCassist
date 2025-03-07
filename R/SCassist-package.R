#' @keywords internal
#' @title SCassist: An AI-Based Workflow Assistant for Single Cell Analysis
#'
#' @description SCassist is a comprehensive toolkit designed to empower single-cell RNA sequencing (scRNA-seq) researchers with intelligent insights and recommendations. The package leverages the power of large language models (LLMs) to streamline the analysis workflow and guide decision-making, ultimately leading to more efficient and interpretable results. 
#' 
#' @details
#' SCassist offers functions for:
#' 
#' Recommending filtering cutoffs.
#' Suggesting appropriate normalization methods. 
#' Summarizing enriched concepts among variable features.
#' Recommending optimal number of PCs to use.
#' Analyzing PCs to glean insights. 
#' Recommending k.param and resolution range. 
#' Predicting cell types.
#' Performing, analyzing, integrating and summarizing KEGG and GO enrichment.
#' 
#' @author Vijay Nagarajan PhD, NEI/NIH
#' @rdname SCassist
#' @seealso  \code{\link{SCassist_analyze_quality}}, 
#'           \code{\link{SCassist_recommend_normalization}}, 
#'           \code{\link{SCassist_analyze_variable_features}},
#'           \code{\link{SCassist_analyze_enrichment}},
#'           \code{\link{SCassist_analyze_pcs}},
#'           \code{\link{SCassist_recommend_k}},
#'           \code{\link{SCassist_recommend_normalization}},
#'           \code{\link{SCassist_recommend_pcs}},
#'           \code{\link{SCassist_recommend_res}},
#'           \code{\link{SCassist_analyze_and_annotate}}
#' @import rollama
#' @importFrom stats sd
## usethis namespace: start
## usethis namespace: end
#' @name SCassist
NULL