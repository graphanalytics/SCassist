# SCassist: An AI-Powered Workflow Assistant for Single-Cell Analysis

**SCassist** utilizes a combination of statistical calculations and LLM-based insights to guide users through the complex process of single-cell RNA-seq data analysis. The package aims to provide recommendations, annotations, and interpretations, leading to efficient and insightful results.

**Features:**

* Recommendations for quality control and filtering
* Recommendations for the most appropriate Normalization method
* Analysis of variable features to identify processes driving the cell-to-cell variation
* Analysis of top PCs to glean insights in to processes driving the system
* Recommendations for optimal number of PCs to use for downstream analysis
* Recommendations for suitable range of resolution values for clustering
* Marker gene analysis and cell type prediction with detailed reasoning
* KEGG pathway and GO enrichment analysis and integration, providing deeper insights in to system understanding

**Benefits:**

* **Automated Recommendations:** Receive tailored recommendations for key parameters and analysis choices based on your specific data characteristics.
* **LLM-Powered Insights:** Leverage the power of LLMs to interpret complex data, uncover hidden patterns, and generate insightful summaries and explanations.
* **Intuitive Interface:** Easy-to-use functions with clear documentation and examples, making advanced analysis accessible to researchers of all levels.
* **Confidential:** Option to use a local LLM server to keep your data and analysis confidential.
* **Cost effective:** If taking advantage of google models, use pay-as-you-go low cost API call options.

**Installation:**

```R
# Install the devtools package if you don't have it
install.packages("devtools")

# Install SCassist from GitHub
devtools::install_github("NIH-NEI/SCassist")
```

**Example Usage:**

**Download example data:** [NK, CD4+ and CD8+ T cells from LCMV infected Ifng - CTCF binding site mutant mice - GSM6625298_scRNA_LCMV_Day4_CD4_CD8_NK_WT_filtered_feature_bc_matrix.h5](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6625299&format=file&file=GSM6625299%5FscRNA%5FLCMV%5FDay4%5FCD4%5FCD8%5FNK%5FKO%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5)


```R
# Load the SCassist and Seurat packages
library(SCassist)
library(Seurat)

# Load the downloaded example file
KO <- Read10X_h5("GSM6625298_scRNA_LCMV_Day4_CD4_CD8_NK_WT_filtered_feature_bc_matrix.h5", use.names = T)

# Create seurat object
KO <- CreateSeuratObject(counts = KO[["Gene Expression"]], names.field = 2,names.delim = "\\-")

# Set api_key_file variable
api_key_file = "api_key_from_google.txt"

# Recommend quality control filters using Gemini (online)
qc_recommendations <- SCassist_analyze_quality("KO", llm_server="google", api_key_file = api_key_file)

# Recommend quality control filters using Llama3 (local)
qc_recommendations <- SCassist_analyze_quality("KO", llm_server="ollama")

# ...and many more functions!
```

**Documentation:**

Detailed documentation for each function, including parameters, usage, and expected outputs, is available through the `?` help function in R.

**License:**

The license for this package can be found in the `LICENSE` file within the package directory.
