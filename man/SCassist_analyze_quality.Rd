% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_quality.R
\name{SCassist_analyze_quality}
\alias{SCassist_analyze_quality}
\title{Analyze Seurat Object Quality and Recommend Filtering Cutoffs}
\usage{
SCassist_analyze_quality(llm_server="google",
                         seurat_object_name,
                         percent_mt = NULL, 
                         percent_ribo = NULL, 
                         percent_hb = NULL,
                         temperature = 0,
                         max_output_tokens = 10048,
                         model_G = "gemini-1.5-flash-latest",
                         model_O = "llama3",
                         api_key_file = "api_keys.txt",
                         model_params = list(seed = 42, temperature = 0, num_gpu = 0))
}
\arguments{
\item{llm_server}{The LLM server to use. Options are "google" or "ollama". Default is "google".}

\item{seurat_object_name}{Character string representing the name of the
Seurat object in the current environment.}

\item{percent_mt}{Numeric value representing the percentage of mitochondrial
genes in each cell. If provided, the function will analyze this metric for
quality control.}

\item{percent_ribo}{Numeric value representing the percentage of ribosomal
genes in each cell. If provided, the function will analyze this metric for
quality control.}

\item{percent_hb}{Numeric value representing the percentage of hemoglobin
genes in each cell. If provided, the function will analyze this metric for
quality control.}

\item{temperature}{Numeric value between 0 and 1, controlling the
"creativity" of the Gemini model's response. Higher values lead to more
diverse and potentially unexpected outputs. Default is 0.}

\item{max_output_tokens}{Integer specifying the maximum number of tokens
the Gemini model can generate in its response. Default is 10048.}

\item{model_G}{Character string specifying the Gemini model to use for
analysis. Default is "gemini-1.5-flash-latest".}

\item{model_O}{Character string specifying the Ollama model to use for
analysis. Default is "llama3".}

\item{api_key_file}{Character string specifying the path to 
a file containing the API key for accessing the Gemini model.}

\item{model_params}{A list of parameters to be passed to the `ollama::query` function.
This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.}
}
\value{
Character string containing the Gemini model's recommendations for
  filtering cutoffs based on the provided data.
}
\description{
This function analyzes a Seurat object and recommends 
filtering cutoffs for quality control based on various metrics, 
including nCount_RNA, nFeature_RNA and user-specified quality 
metrics (e.g., percent mitochondrial, ribosomal or hemoglobin genes). 
It leverages the Gemini language model to provide insightful 
recommendations based on the data distribution and potential
outliers.
}
\examples{
\dontrun{
# Assuming you have a Seurat object named 'seurat_obj'
SCassist_analyze_quality(seurat_object_name = "seurat_obj",
                         percent_mt = "percent.mt",
                         percent_ribo = "percent.ribo",
                         percent_hb = "percent.hb",
                         api_key_file = "my_api_key.txt")
}

}
