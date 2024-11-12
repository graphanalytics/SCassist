#' @title Recommend Number of Principal Components (PCs) Using an LLM
#'
#' @description This function analyzes the variance explained by each principal component
#' (PC) from a PCA analysis using a large language model (LLM) to identify
#' the "elbow point" in the scree plot. Based on this elbow point, the LLM
#' recommends the number of PCs to use for downstream analysis, such as
#' finding neighbors or running UMAP.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#' 
#' @details This function was written with assistance from Google's Gemini and Meta's Llama3.
#'
#' @param llm_server The LLM server to use. Options are "google" or "ollama". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name The name of the Seurat object containing the
#'  single-cell RNA-seq data. The object should be accessible in the current
#'  environment and should have PCA already run (e.g., using `RunPCA`).
#' @param experimental_design (Optional) A character string describing the experimental
#'  design. This information helps the LLM contextualize the analysis. If not provided,
#'  the prompt will not include the experimental design information.
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
#' @return A numeric value representing the recommended number of PCs to use
#'  for downstream analysis.
#'
#' @usage SCassist_recommend_pcs(llm_server="google",
#'  seurat_object_name, 
#'  experimental_design = NULL,
#'  temperature = 0,
#'  max_output_tokens = 10048,
#'  model_G = "gemini-1.5-flash-latest",
#'  model_O = "llama3",
#'  api_key_file = "api_keys.txt",
#'  model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#'  
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj' with PCA run
#' SCassist_recommend_pcs(seurat_object_name = "seurat_obj",
#'                        experimental_design = "Time course experiment",
#'                        api_key_file = "my_api_key.txt")
#' }
#' @import rollama
#' @keywords single-cell, PCA, dimensionality reduction, LLM, recommendation
#' @export

SCassist_recommend_pcs <- function(llm_server="google",
                                     seurat_object_name,
                                     experimental_design = NULL,
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_recommend_pcs_G(seurat_object_name,
                                      experimental_design = experimental_design,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_recommend_pcs_L(seurat_object_name, 
                                      experimental_design = experimental_design,
                                      model_params = model_params, 
                                      model = model_O ))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_recommend_pcs_G <- function(seurat_object_name, experimental_design = NULL,
                                   temperature = 0,
                                   max_output_tokens = 10048,
                                   model = "gemini-1.5-flash-latest",
                                   api_key_file = "api_keys.txt") {
  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # 2. Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, 
           "' not found in environment.", call. = FALSE)
    }
  )
  
  # 3. Check if PCA is available
  if (is.null(seurat_object@reductions$pca)) {
    stop("Please run PCA on the Seurat object first (e.g., using `RunPCA`).")
  }
  
  # 4. Prepare Data for LLM
  variance_explained <- seurat_object@reductions$pca@stdev^2
  variance_explained <- variance_explained / sum(variance_explained) * 100
  variance_text <- paste0("PC", 1:length(variance_explained), ": ", 
                          round(variance_explained, 2), "%", collapse = ", ")
  
  # 5. Craft a Prompt
  prompt_start <- "I have a single-cell experiment where I performed principal component analysis (PCA). The variance explained by each PC is: "
  
  if (!is.null(experimental_design)) {
    prompt_start <- paste0(
      "I have a single-cell experiment where I performed principal component analysis (PCA) on ", 
      experimental_design, ". The variance explained by each PC is:

"
    )
  }
  
  prompt <- paste0(
    prompt_start,
    variance_text, 
    
    "Based on this information, determine the optimal number of PCs to use for downstream analysis, such as finding neighbors or running UMAP. Explain your reasoning and consider the following factors:

* **Elbow point:** Is there a clear 'elbow' in the scree plot?
* **Variance explained:**  What percentage of the total variance is captured by the chosen number of PCs?
* **Balance between complexity and interpretability:** A higher number of PCs might capture more subtle variation but make the analysis more complex.

Provide your recommendation as a single number (e.g., 5) and a concise explanation.
")
  
  # 6. Send the request to the Gemini API
  response1 <- POST(
    url = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model_query),
    query = list(key = api_key), # Access api_key from api_keys.R
    content_type_json(),
    encode = "json",
    body = list(
      contents = list(
        parts = list(
          list(text = prompt)
        )),
      generationConfig = list(
        temperature = temperature,
        maxOutputTokens = max_output_tokens,
        seed = 123456  # For reproducible runs
      )
    )
  )
  
  # 7. Extract the Gemini model's response
  candidates <- content(response1)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  
  # 8. Print and return the LLM's response
  cat(outputs[["text"]])
  return(outputs[["text"]])
}

SCassist_recommend_pcs_L <- function(seurat_object_name, experimental_design = NULL,
                                     model_params = list(seed = 42, temperature = 0),
                                     model = "llama3") {
  # Get the Seurat object
  seurat_object <- tryCatch(
    {
      get(seurat_object_name)
    },
    error = function(e) {
      stop("Error: Seurat object '", seurat_object_name, 
           "' not found in environment.", call. = FALSE)
    }
  )
  
  # Check if PCA is available
  if (is.null(seurat_object@reductions$pca)) {
    stop("Please run PCA on the Seurat object first (e.g., using `RunPCA`).")
  }
  
  # Prepare Data for LLM
  variance_explained <- seurat_object@reductions$pca@stdev^2
  variance_explained <- variance_explained / sum(variance_explained) * 100
  variance_text <- paste0("PC", 1:length(variance_explained), ": ", 
                          round(variance_explained, 2), "%", collapse = ", ")
  
  # Craft a Prompt
  prompt_start <- "I have a single-cell experiment where I performed principal component analysis (PCA). The variance explained by each PC is:

"
  
  if (!is.null(experimental_design)) {
    prompt_start <- paste0(
      "I have a single-cell experiment where I performed principal component analysis (PCA) on ", 
      experimental_design, ". The variance explained by each PC is:

"
    )
  }
  
  prompt <- paste0(
    prompt_start,
    variance_text, 
    
    "Based on this information, determine the optimal number of PCs to use for downstream analysis, such as finding neighbors or running UMAP. Explain your reasoning and consider the following factors:

* **Elbow point:** Is there a clear 'elbow' in the scree plot?
* **Variance explained:**  What percentage of the total variance is captured by the chosen number of PCs?
* **Balance between complexity and interpretability:** A higher number of PCs might capture more subtle variation but make the analysis more complex.

Provide your recommendation as a single number (e.g., 5) and a concise explanation.
")
  
  # Query the LLM
  response <- tryCatch(
    {
      suppressMessages(
        rollama::query(
          prompt,
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
  if (is.null(response$message$content)) {
    stop("Error: The LLM returned an invalid response. Please check the LLM model and parameters.")
  }
  
  # Extract the recommended number of PCs
  recommended_pcs_text <- gsub(".*?([0-9]+).*", "\\1", response$message$content)
  recommended_pcs <- as.numeric(recommended_pcs_text)
  
  # Handle invalid response
  if (is.na(recommended_pcs) || recommended_pcs <= 0) {
    warning("The LLM could not provide a valid recommendation. Please check the LLM model and parameters.")
    recommended_pcs <- 0
  }
  
  # Print and return the LLM's response
  cat(response$message$content)
  return(recommended_pcs)
}