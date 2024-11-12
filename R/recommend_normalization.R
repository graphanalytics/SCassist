#' @title Recommend Normalization Method for Single-Cell RNA-seq Data
#'
#' @description This function analyzes key characteristics of a single-cell RNA-seq dataset
#' using a large language model (LLM) to recommend the most appropriate
#' normalization method from Seurat's available options. The LLM considers factors
#' like the number of cells, gene expression distribution, and library size
#' variation to provide a reasoned recommendation.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#'
#' @details This function was written with assistance from 
#' Google's Gemini and Meta's Llama3. 
#' 
#' @usage 
#' SCassist_recommend_normalization(llm_server="google",
#'  seurat_object_name, 
#'  temperature = 0,
#'  max_output_tokens = 10048,
#'  model_G = "gemini-1.5-flash-latest",
#'  model_O = "llama3",
#'  api_key_file = "api_keys.txt",
#'  model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#'  
#' @param llm_server The LLM server to use. Options are "google" or "ollama". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name The name of the Seurat object containing the
#'  single-cell RNA-seq data. The object should be accessible in the current
#'  environment.
#' @param temperature Numeric value between 0 and 1, controlling the
#'   "creativity" of the Gemini model's response. Higher values lead to more
#'   diverse and potentially unexpected outputs. Default is 0.
#' @param max_output_tokens Integer specifying the maximum number of tokens
#'   the Gemini model can generate in its response. Default is 10048.
#' @param model_G Character string specifying the Gemini model to use for
#'  analysis. Default is "gemini-1.5-flash-latest".
#' @param model_O Character string specifying the Ollama model to use for
#'  analysis. Default is "llama3".
#' @param api_key_file Character string specifying the path to a file containing
#'   the API key for accessing the Gemini model.
#'
#' @return A character string containing the LLM's recommendation for the
#'  best normalization method, along with a justification for the choice and 
#'  discussion of potential alternatives.
#'
#' @import rollama
#' @importFrom stats sd
#' @importFrom Seurat DefaultAssay
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj'
#' SCassist_recommend_normalization(seurat_object_name = "seurat_obj",
#'                                  api_key_file = "my_api_key.txt")
#' }
#'
#' @export

SCassist_recommend_normalization <- function(llm_server="google",
                                     seurat_object_name,
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_recommend_normalization_G(seurat_object_name,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_recommend_normalization_L(seurat_object_name, 
                                      model_params = model_params, 
                                      model = model_O ))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_recommend_normalization_G <- function(seurat_object_name, 
                                               temperature = 0,
                                               max_output_tokens = 10048,
                                               model = "gemini-1.5-flash-latest",
                                               api_key_file = "api_keys.txt") {
  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # 2. Ensure Seurat is available (using requireNamespace)
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for this function. Please install it.", call. = FALSE)
  }
  
  # 3. Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # 4. Dummy function call to satisfy the dependency check
  Seurat::DefaultAssay(seurat_object)
  
  # 5. Extract relevant data
  num_cells <- ncol(seurat_object)
  gene_expr_mean <- mean(seurat_object@assays$RNA$counts@x)  # Access count values
  gene_expr_sd <- sd(seurat_object@assays$RNA$counts@x)
  library_sizes <- Matrix::colSums(seurat_object@assays$RNA$counts)
  library_size_cv <- sd(library_sizes) / mean(library_sizes) 
  
  # 6. Get available normalization methods from Seurat
  seurat_normalization_methods <- c("LogNormalize", "CLR", "RC", "SCTransform")
  
  # 7. Prepare input for LLM
  input_text <- paste0("Consider the specific characteristics of the single-cell RNA-seq dataset, such as the number of cells, the distribution of gene expression values, and the observed library size variations. Recommend the most appropriate normalization method for this dataset from the following options available in Seurat: ",
                       paste(seurat_normalization_methods, collapse = ", "), ". Explain why this method is preferred and discuss any potential alternatives.) \n\n",
                       "The single-cell RNA-seq dataset contains ", num_cells, " cells.  The mean gene expression is ", gene_expr_mean, " with a standard deviation of ", gene_expr_sd, ". The coefficient of variation for library sizes is ", library_size_cv, ".")
  
  # 8. Send the request to the Gemini API
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
  
  # 9. Extract the Gemini model's response
  candidates <- content(response1)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  
  # 10. Print and return the LLM's response
  cat(outputs[["text"]])
  return(outputs[["text"]])
}

SCassist_recommend_normalization_L <- function(seurat_object_name, model_params = list(seed = 42, temperature = 0), 
                                               model = "llama3") {
  
  # Ensure Seurat is available (using requireNamespace)
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for this function. Please install it.", call. = FALSE)
  }
  
  # Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, "' not found in environment.", call. = FALSE)
    }
  )
  
  # Dummy function call to satisfy the dependency check
  Seurat::DefaultAssay(seurat_object)
  
  # Extract relevant data
  num_cells <- ncol(seurat_object)
  gene_expr_mean <- mean(seurat_object@assays$RNA$counts@x)  # Access count values
  gene_expr_sd <- sd(seurat_object@assays$RNA$counts@x)
  library_sizes <- Matrix::colSums(seurat_object@assays$RNA$counts)
  library_size_cv <- sd(library_sizes) / mean(library_sizes) 
  
  # Get available normalization methods from Seurat
  seurat_normalization_methods <- c("LogNormalize", "CLR", "RC", "SCTransform")
  
  # Prepare input for LLM
  input_text <- paste0("Consider the specific characteristics of the single-cell RNA-seq dataset, such as the number of cells, the distribution of gene expression values, and the observed library size variations. Recommend the most appropriate normalization method for this dataset from the following options available in Seurat: ",
                       paste(seurat_normalization_methods, collapse = ", "), ". Explain why this method is preferred and discuss any potential alternatives.) \n\n",
                       "The single-cell RNA-seq dataset contains ", num_cells, " cells.  The mean gene expression is ", gene_expr_mean, " with a standard deviation of ", gene_expr_sd, ". The coefficient of variation for library sizes is ", library_size_cv, ".")
  
  # Query the LLM
  response <- tryCatch(
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
           "Please check your ollama server connection and/or the model.", call. = FALSE)
    }
  )
  
  # Check if LLM response is valid
  if (is.null(response$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
  }
  
  # Print and return the LLM's response
  cat(response$message$content) 
  return(response$message$content)
}