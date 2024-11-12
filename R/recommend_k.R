#' @title Recommend `k.param` Value for Seurat's `FindNeighbors` Function
#'
#' @description This function analyzes key characteristics of a single-cell RNA-seq dataset
#' and the desired clustering goals using a large language model (LLM) to
#' recommend a range of potential `k.param` values for Seurat's `FindNeighbors`
#' function. The LLM considers factors like the number of cells, the number of
#' principal components (PCs) used for dimensionality reduction, and the
#' experimental design and clustering goals to provide a reasoned suggestion.
#'
#' @author Vijay Nagarajan, PhD, NEI/NIH
#'
#' @details This function was written with assistance from Google's Gemini and Meta's Llama3.
#'
#' @usage SCassist_recommend_k(llm_server="google", 
#'                  seurat_object_name, num_pcs, 
#'                  experimental_design = "",
#'                  temperature = 0,
#'                  max_output_tokens = 10048,
#'                  model_G = "gemini-1.5-flash-latest",
#'                  model_O = "llama3",
#'                  api_key_file = "api_keys.txt",
#'                  model_params = list(seed = 42, temperature = 0, num_gpu = 0))
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj'
#' SCassist_recommend_k(seurat_object_name = "seurat_obj",
#'                      num_pcs = 20,
#'                      experimental_design = "Single-cell RNA-seq data from a mouse brain",
#'                      api_key_file = "my_api_key.txt")
#' }
#' @param llm_server The LLM server to use. Options are "google" or "ollama". Default is "google".
#' @param model_params A list of parameters to be passed to the `ollama::query` function.
#'   This allows customization of the Llama model's behavior. Default is `list(seed = 42, temperature = 0, num_gpu = 0)`.
#' @param seurat_object_name The name of the Seurat object containing the
#'  single-cell RNA-seq data. The object should be accessible in the current
#'  environment.
#' @param num_pcs The number of principal components (PCs) used for
#'  dimensionality reduction.
#' @param experimental_design (Optional) A character string describing the
#'  experimental design. This information helps the LLM contextualize the
#'  analysis. It should include information about the sample source,
#'  experimental conditions, research question, and clustering goals. If not
#'  provided, the LLM will analyze the data without specific experimental
#'  context.
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
#' @return A character string containing the LLM's recommendation for a range
#'  of `k.param` values, along with a justification for the choice. The
#'  output will be structured as:
#'  - "Recommended K:" followed by a range of whole number values.
#'  - "Reasoning:" followed by a paragraph explaining the rationale behind
#'  the recommendation.
#'
#' @import rollama
#' @export

SCassist_recommend_k <- function(llm_server="google",
                                     seurat_object_name,
                                     num_pcs,
                                     experimental_design = "",
                                     temperature = 0,
                                     max_output_tokens = 10048,
                                     model_G = "gemini-1.5-flash-latest",
                                     model_O = "llama3",
                                     api_key_file = "api_keys.txt",
                                     model_params = list(seed = 42, temperature = 0, num_gpu = 0)
) {
  
  if (llm_server == "google") {
    return(SCassist_recommend_k_G(seurat_object_name,
                                      num_pcs = num_pcs,
                                      experimental_design = experimental_design,
                                      temperature = temperature,
                                      max_output_tokens = max_output_tokens,
                                      model = model_G,
                                      api_key_file = api_key_file ))
  } else if (llm_server == "ollama") {
    return(SCassist_recommend_k_L(seurat_object_name, 
                                      num_pcs = num_pcs,
                                      experimental_design = experimental_design,
                                      model_params = model_params, 
                                      model = model_O ))
  } else {
    stop("Invalid llm_server option. Please specify 'google' or 'ollama'.")
  }
}

SCassist_recommend_k_G <- function(seurat_object_name, num_pcs,
                                 experimental_design = "", 
                                 temperature = 0,
                                 max_output_tokens = 10048,
                                 model = "gemini-1.5-flash-latest",
                                 api_key_file = "api_keys.txt") {
  
  # 1. Set up the Gemini API request
  model_query <- paste0(model, ":generateContent")
  
  # Read the API key from the specified file
  api_key <- readLines(api_key_file)
  
  # 2. Check if num_pcs is provided
  if (is.null(num_pcs)) {
    stop("Error: Please provide the number of principal components (num_pcs) used for dimensionality reduction.", call. = FALSE)
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
  
  # 4. Extract relevant data
  num_cells <- ncol(seurat_object)
  
  # 5. Set research question and clustering goal
  research_question <- "What are the different cell populations present in the sample?"
  clustering_goal <- "identify biologically meaningful clusters representing the diverse cell types in the sample"
  
  # 6. Construct the prompt for the LLM
  prompt <- paste0("I'm analyzing single-cell RNA sequencing data from \n", experimental_design, "\n",
                   "The dataset contains approximately ", num_cells, " cells \n",
                   "I've determined that using ", num_pcs, " PCs (`dims` = ", num_pcs, ") is suitable for my data. ",
                   "I'm interested in identifying distinct cell populations.",
                   "My goal is to", clustering_goal, ". \n",
                   "Can you suggest a range of potential `k.param` values, in whole number, to explore for `FindNeighbors()` based on this information? Provide the output as, Recommended K: and a two short reasoning paragraph under Reasoning: ")
  
  # 7. Send the request to the Gemini API
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
  
  # 8. Extract the Gemini model's response
  candidates <- content(response1)$candidates
  outputs <- unlist(lapply(candidates, function(candidate) candidate$content$parts))
  
  # 9. Print and return the LLM's response
  cat(outputs[["text"]])
  return(outputs[["text"]])
}

SCassist_recommend_k_L <- function(seurat_object_name, num_pcs,
                                   experimental_design = "", model_params = list(seed = 42, temperature = 0), model = "llama3") {
  
  # Check if num_pcs is provided
  if (is.null(num_pcs)) {
    stop("Error: Please provide the number of principal components (num_pcs) used for dimensionality reduction.", call. = FALSE)
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
  
  # Extract relevant data
  num_cells <- ncol(seurat_object)
  
  # Set research question and clustering goal
  research_question <- "What are the different cell populations present in the sample?"
  clustering_goal <- "identify biologically meaningful clusters representing the diverse cell types in the sample"
  
  # Construct the prompt for the LLM
  prompt <- paste0("I'm analyzing single-cell RNA sequencing data from \n", experimental_design, "\n",
                   "The dataset contains approximately ", num_cells, " cells \n",
                   "I've determined that using ", num_pcs, " PCs (`dims` = ", num_pcs, ") is suitable for my data. ",
                   "I'm interested in identifying distinct cell populations.",
                   "My goal is to", clustering_goal, ". \n",
                   "Can you suggest a range of potential `k.param` values, in whole number, to explore for `FindNeighbors()` based on this information? Provide the output as, Recommended K: and a two short reasoning paragraph under Reasoning: ")
  
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
  
  # Print the LLM's response
  cat(response$message$content)
  
  # Return the LLM's response (for potential use elsewhere)
  return(response$message$content)
}