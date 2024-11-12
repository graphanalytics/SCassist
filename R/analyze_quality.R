#' @title Analyze Seurat Object Quality and Recommend Filtering Cutoffs
#'
#' @description This function analyzes a Seurat object and recommends 
#' filtering cutoffs for quality control based on various metrics, 
#' including nCount_RNA, nFeature_RNA and user-specified quality 
#' metrics (e.g., percent mitochondrial, ribosomal or hemoglobin genes). 
#' It leverages the Gemini language model to provide insightful 
#' recommendations based on the data distribution and potential
#' outliers.
#'
#' @param llm_server The LLM server to use. Options are "google" or "ollama". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name Character string representing the name of the
#'   Seurat object in the current environment.
#' @param percent_mt Numeric value representing the percentage of mitochondrial
#'   genes in each cell. If provided, the function will analyze this metric for
#'   quality control.
#' @param percent_ribo Numeric value representing the percentage of ribosomal
#'   genes in each cell. If provided, the function will analyze this metric for
#'   quality control.
#' @param percent_hb Numeric value representing the percentage of hemoglobin
#'   genes in each cell. If provided, the function will analyze this metric for
#'   quality control.
#' @param temperature Numeric value between 0 and 1, controlling the
#'   "creativity" of the Gemini model's response. Higher values lead to more
#'   diverse and potentially unexpected outputs. Default is 0.
#' @param max_output_tokens Integer specifying the maximum number of tokens
#'   the Gemini model can generate in its response. Default is 10048.
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#' @param api_key_file Character string specifying the path to 
#'  a file containing the API key for accessing the Gemini model.
#'
#' @return Character string containing the Gemini model's recommendations for
#'   filtering cutoffs based on the provided data.
#'
#' @usage
#' SCassist_analyze_quality(llm_server="google",
#'                          seurat_object_name,
#'                          percent_mt = NULL, 
#'                          percent_ribo = NULL, 
#'                          percent_hb = NULL,
#'                          temperature = 0,
#'                          max_output_tokens = 10048,
#'                          model_G = "gemini-1.5-flash-latest",
#'                          model_O = "llama3",
#'                          api_key_file = "api_keys.txt",
#'                          model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj'
#' SCassist_analyze_quality(seurat_object_name = "seurat_obj",
#'                          percent_mt = "percent.mt",
#'                          percent_ribo = "percent.ribo",
#'                          percent_hb = "percent.hb",
#'                          api_key_file = "my_api_key.txt")
#' }
#'
#' @importFrom httr POST content content_type_json
#' @importFrom stats quantile
#'
#' @export

SCassist_analyze_quality <- function(llm_server="google",
                                     seurat_object_name,
                                     percent_mt = NULL, 
                                     percent_ribo = NULL, 
                                     percent_hb = NULL,
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_analyze_quality_G(seurat_object_name,
                                      percent_mt = percent_mt, 
                                      percent_ribo = percent_ribo, 
                                      percent_hb = percent_hb,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_analyze_quality_L(seurat_object_name, 
                                      model_params = model_params, 
                                      model = model_O, 
                                      percent_mt = percent_mt, 
                                      percent_ribo = percent_ribo, 
                                      percent_hb = percent_hb ))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_analyze_quality_G <- function(seurat_object_name,
                                       percent_mt = NULL, 
                                       percent_ribo = NULL, 
                                       percent_hb = NULL,
                                       temperature = 0,
                                       max_output_tokens = 10048,
                                       model = "gemini-1.5-flash-latest",
                                       api_key_file = "api_keys.txt") {
  
  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # 2. Retrieve the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # 3. Identify available quality metrics
  quality_metrics <- if (!is.null(percent_mt) || !is.null(percent_ribo) || !is.null(percent_hb)) {
    intersect(c(percent_mt, percent_ribo, percent_hb), colnames(seurat_object@meta.data))
  } else {
    NULL
  }
  
  # 4. Calculate summary and quantile data for nCount_RNA and nFeature_RNA
  nCount_summary <- summary(seurat_object$nCount_RNA)
  nCount_quantiles <- quantile(seurat_object$nCount_RNA, probs = c(0.05, 0.95))
  
  nFeature_summary <- summary(seurat_object$nFeature_RNA)
  nFeature_quantiles <- quantile(seurat_object$nFeature_RNA, probs = c(0.05, 0.95))
  
  # 5. Build input text for the Gemini model
  input_text <- paste0(
    "I have summary statistics and quantile data for nCount_RNA, nFeature_RNA, ",
    paste(quality_metrics, collapse = ", "), ", from a single cell experiment.
    I want to determine refined cutoff values that are more sensitive to the tails of each distribution, by combining both the summary statistics and the quantile information. The goal is to filter out poor quality cells.

  **nCount_RNA**

  **Summary Statistics:**
  * Minimum: ", nCount_summary["Min."], "\n",
    "* 1st Quartile: ", nCount_summary["1st Qu."], "\n",
    "* Median: ", nCount_summary["Median"], "\n",
    "* Maximum: ", nCount_summary["Max."], "\n\n",
    "**Quantile Data:**
  * 5th Percentile: ", nCount_quantiles["5%"], "\n",
    "* 95th Percentile: ", nCount_quantiles["95%"], "\n\n",
    
    "**nFeature_RNA**

  **Summary Statistics:**
  * Minimum: ", nFeature_summary["Min."], "\n",
    "* 1st Quartile: ", nFeature_summary["1st Qu."], "\n",
    "* Median: ", nFeature_summary["Median"], "\n",
    "* Maximum: ", nFeature_summary["Max."], "\n\n",
    "**Quantile Data:**
  * 5th Percentile: ", nFeature_quantiles["5%"], "\n",
    "* 95th Percentile: ", nFeature_quantiles["95%"], "\n\n"
  )
  
  # 6. Calculate and print summary statistics for quality metrics
  all_output <- capture.output({
    for (metric in quality_metrics) {
      # Extract the data as a vector
      data_vector <- seurat_object@meta.data[[metric]]
      
      # Calculate summary statistics (using 'summary' for non-zero values)
      nonzero_values <- data_vector[data_vector != 0]
      summary_data <- summary(nonzero_values)
      
      # Calculate 95th percentile
      quantile_data <- quantile(data_vector, probs = 0.95)
      
      # Capture the printed output of 'summary_data' as a string
      summary_output <- capture.output(summary_data)
      
      # Extract the header and values
      header <- summary_output[1]
      values <- paste(summary_output[-1], collapse = " ")
      
      # Combine the header and values into a single string with proper formatting
      summary_output <- paste0(header, "\n", values)
      
      # Print the summary statistics (formatted)
      cat(paste0("Summary statistics for ", metric, ":\n", summary_output, "\n\n"))
      
      # Print the 95th percentile
      cat(paste0("95th percentile data for ", metric, ": ", quantile_data, "\n\n"))
    }
  })
  
  # 7. Combine all output into a single string
  all_output <- paste(all_output, collapse = "\n")
  
  # 8. Combine the strings, adding a newline between them
  combined_text <- paste(input_text, all_output, sep = "\n\n")
  
  # 9. Get the total number of cells
  nCells <- length(seurat_object@meta.data$orig.ident)
  
  # 10. Complete the input text for the Gemini model
  input_text <- paste0(
    combined_text,
    "** Total number of cells : **", nCells, "\n",
    "**Please provide the refined lower and upper cutoff values for nCount_RNA and nFeature_RNA, ",
    ifelse(!is.null(quality_metrics), paste0("and the upper cutoff values for ", paste(quality_metrics, collapse = ", "), ", "), ""),
    "calculated based on the provided data and taking into account the potential presence of tails in the distribution. Dont generate the cuffoff ONLY based on the percentiles, instead combine the percentiles with a logic that we need to remove outliers, but still keep good data. Make sure to state that the researcher should test a range of values around your recommendation**. Do not provide any IMPORTANT NOTES. Start the response saying that, based on the your data summary, below are my recommendations for the quality filtering of the data"
  )
  
  # 11. Send the request to the Gemini API
  response1 <- POST(
    url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
    query = list(key = api_key), # Access api_key from api_keys.R
    content_type_json(),
    encode = "json",
    body = list(
      contents = list(
        parts = list(
          list(text = input_text)
        )),
      generationConfig = list(
        temperature = temperature,
        maxOutputTokens = max_output_tokens,
        seed = 123456  # For reproducible runs
      )
    )
  )
  
  # 12. Extract the Gemini model's response
  candidates <- content(response1)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  
  # 13. Print and return the LLM's response
  cat(outputs[["text"]])
  return(outputs[["text"]])
}


SCassist_analyze_quality_L <- function(seurat_object_name, 
                                       percent_mt = NULL, 
                                       percent_ribo = NULL, 
                                       percent_hb = NULL,
                                       model_params = list(seed = 42, temperature = 0), 
                                       model = "llama3") {
  
  # Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # Identify available quality metrics
  quality_metrics <- if (!is.null(percent_mt) || !is.null(percent_ribo) || !is.null(percent_hb)) {
    intersect(c(percent_mt, percent_ribo, percent_hb), colnames(seurat_object@meta.data))
  } else {
    NULL
  }
  
  
  # Calculate summary and quantile data
  nCount_summary <- summary(seurat_object$nCount_RNA)
  nCount_quantiles <- quantile(seurat_object$nCount_RNA, probs = c(0.05, 0.95))
  
  nFeature_summary <- summary(seurat_object$nFeature_RNA)
  nFeature_quantiles <- quantile(seurat_object$nFeature_RNA, probs = c(0.05, 0.95))
  
  # Build input text with only available metrics
  input_text <- paste0(
    "I have summary statistics and quantile data for nCount_RNA, nFeature_RNA, ", 
    paste(quality_metrics, collapse = ", "), ", from a single cell experiment. 
    I want to determine refined cutoff values that are more sensitive to the tails 
    of each distribution, by combining both the summary statistics and the quantile     information. The goal is to filter out poor quality cells.
  
  **nCount_RNA**

  **Summary Statistics:**
  * Minimum: ", nCount_summary["Min."], "\n",
    "* 1st Quartile: ", nCount_summary["1st Qu."], "\n",
    "* Median: ", nCount_summary["Median"], "\n",
    "* Maximum: ", nCount_summary["Max."], "\n\n",
    "**Quantile Data:**
  * 5th Percentile: ", nCount_quantiles["5%"], "\n",
    "* 95th Percentile: ", nCount_quantiles["95%"], "\n\n",
    
    "**nFeature_RNA**

  **Summary Statistics:**
  * Minimum: ", nFeature_summary["Min."], "\n",
    "* 1st Quartile: ", nFeature_summary["1st Qu."], "\n",
    "* Median: ", nFeature_summary["Median"], "\n",
    "* Maximum: ", nFeature_summary["Max."], "\n\n",
    "**Quantile Data:**
  * 5th Percentile: ", nFeature_quantiles["5%"], "\n",
    "* 95th Percentile: ", nFeature_quantiles["95%"], "\n\n"
  )
  
  # Add data for each available quality metric
  for (metric in quality_metrics) {
    # Extract the vector from the Seurat object's metadata
    metric_data <- seurat_object@meta.data[[metric]]
    
    # Handle cases where the metric is a single vector
    if (is.vector(metric_data)) {
      summary_data <- summary(metric_data)
      quantile_data <- quantile(metric_data, probs = c(0.95))
    } else {
      # If the metric is a data frame, you can use the mean or median
      summary_data <- summary(metric_data) # This will give summary statistics
      quantile_data <- quantile(rowMeans(metric_data), probs = c(0.95)) # Or use rowMeans for median 
    }
    
    input_text <- paste0(
      input_text, 
      "**", metric, "**

  **Summary Statistics:**
  * ", summary_data, "\n",
      "**Quantile Data:**
  * 95th Percentile: ", quantile_data["95%"], "\n\n"
    )
  }
  
  nCells <- length(seurat_object@meta.data$orig.ident)
  
  # Complete the input text
  input_text <- paste0(
    input_text,
    "** Total number of cells : **", nCells, "\n",
    
    "**Please provide the refined lower and upper cutoff values for nCount_RNA and nFeature_RNA, and the upper cutoff values for ", 
    paste(quality_metrics, collapse = ", "), ", calculated based on the provided data and taking into account the potential presence of tails in the distribution. Make sure to state that the researcher should test a range of values around your recommendation**"
  )
  
  # Query the LLM
  response1 <- tryCatch(
    {
      suppressMessages(
        rollama::query(
          input_text,
          model = model,
          screen = FALSE,
          model_params = model_params
        )
      )
    },
    error = function(e) {
      stop("Error: The LLM query encountered an error. ",
           "Please check your ollama server connection and or the model.", call. = FALSE)
    }
  )
  
  # Check if LLM response is valid
  if (is.null(response1$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
  }
  
  # Print and return the LLM's response
  cat(response1$message$content) 
  return(response1$message$content)
}