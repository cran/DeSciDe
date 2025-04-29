## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

## ----setup, include = FALSE---------------------------------------------------
library(DeSciDe)

# Define evaluation flag
eval_flag <- Sys.getenv("NOT_CRAN") == "TRUE"

# Load precomputed results if not in CRAN
if(!eval_flag) {
  load("data/results_default.RData")
  load("data/results_total.RData")
  load("data/threshold_50.RData")
  load("data/threshold_20.RData")
}

## ----error=FALSE--------------------------------------------------------------
# Import genes list and terms list from CSV
genes <- read.csv("genes.csv", header = FALSE)[[1]]
terms <- read.csv("terms.csv", header = FALSE)[[1]]

## -----------------------------------------------------------------------------
genes

## -----------------------------------------------------------------------------
terms

## ----plot_chunk, message=FALSE, warning=FALSE, eval = eval_flag---------------
# results <- descide(genes_list = genes, terms_list = terms)

## ----include = FALSE, eval=!eval_flag-----------------------------------------
# Load precomputed results
results <- results_default

## ----fig.width=8, fig.height=6, echo = FALSE----------------------------------
head(results$summary_results)

## ----fig.width=8, fig.height=6, echo = FALSE----------------------------------
plot_heatmap(results$pubmed_results)

## ----fig.width=8, fig.height=6, echo = FALSE----------------------------------
knitr::include_graphics("data/Network_full.pdf")

## ----fig.width=8, fig.height=6, echo = FALSE----------------------------------
plot_clustering(results$string_results)

## ----fig.width=8, fig.height=6, echo = FALSE----------------------------------
plot_connectivity_precedence(results$summary_results)

## ----message=FALSE, fig.width=8, fig.height=6, eval = eval_flag---------------
# results_total <- descide(genes_list = genes, terms_list = terms, rank_method = "total")

## ----include = FALSE, eval=!eval_flag-----------------------------------------
# Load precomputed results for rank_method = "total"
results_total <- results_total

## -----------------------------------------------------------------------------
head(results$summary_results)

## -----------------------------------------------------------------------------
head(results_total$summary_results)

## ----eval=FALSE, fig.width=8, fig.height=6------------------------------------
# # Change species to mus musculus for STRING search.
# descide(genes_list = genes, terms_list = terms, species = 10090)

## ----eval=FALSE, fig.width=8, fig.height=6------------------------------------
# # Change STRING score threshold to 600.
# descide(genes_list = genes, terms_list = terms, score_threshold = 600)

## ----eval=FALSE, fig.width=8, fig.height=6------------------------------------
# # Change STRING network type to only include physical interactions.
# descide(genes_list = genes, terms_list = terms, network_type = "physical")

## ----full_string_chunk, warning=FALSE, fig.width=8, fig.height=6, eval = eval_flag----
# # Run STRING search and display network with full network.
# full_string <- search_string_db(genes_list = genes, network_type = "full")
# plot_string_network(full_string$string_db, full_string$string_ids)

## ----warning=FALSE, fig.width=8, fig.height=6, echo=FALSE, eval=!eval_flag----
knitr::include_graphics("data/Network_full.pdf")

## ----physical_string_chunk, warning=FALSE, fig.width=8, fig.height=6, eval = eval_flag----
# # Run STRING search and display network with physical network.
# physical_string <- search_string_db(genes_list = genes, network_type = "physical")
# plot_string_network(physical_string$string_db, physical_string$string_ids)

## ----warning=FALSE, fig.width=8, fig.height=6, echo=FALSE, eval=!eval_flag----
knitr::include_graphics("data/Network_physical.pdf")

## ----eval=FALSE, fig.width=8, fig.height=6------------------------------------
# # Command to adjust threshold_percentage for full descide pipeline.
# results <- descide(genes_list = genes, terms_list = terms, threshold_percentage = 50)

## ----fig.width=8, fig.height=6------------------------------------------------
# Calculate and plot threshold of 20%.
threshold_20 <- combine_summary(pubmed_search_results = results$pubmed_results, string_results = results$string_results, threshold_percentage = 20)
plot_connectivity_precedence(combined_summary = threshold_20)

## ----fig.width=8, fig.height=6, include=FALSE, eval=!eval_flag----------------
# Load precomputed results for threshold 20%
threshold_20 <- threshold_20
plot_connectivity_precedence(combined_summary = threshold_20)

## ----fig.width=8, fig.height=6------------------------------------------------
# Calculate and plot threshold of 50%.
threshold_50 <- combine_summary(pubmed_search_results = results$pubmed_results, string_results = results$string_results, threshold_percentage = 50)
plot_connectivity_precedence(combined_summary = threshold_50)

## ----fig.width=8, fig.height=6, include=FALSE, eval=!eval_flag----------------
# Load precomputed results for threshold 50%
threshold_50 <- threshold_50
plot_connectivity_precedence(combined_summary = threshold_50)

## -----------------------------------------------------------------------------
head(threshold_20)
head(threshold_50)

## ----eval=FALSE,fig.width=8, fig.height=6-------------------------------------
# # Code to run DeSciDe and export all plots and tables to desired directory.
# descide(genes_list = genes, terms_list = terms, export = TRUE, file_directory = "your/desired/directory", export_format = "excel")

## ----eval=FALSE---------------------------------------------------------------
# descide(
#   genes_list,
#   terms_list,
#   rank_method = "weighted",
#   species = 9606,
#   network_type = "full",
#   score_threshold = 400,
#   threshold_percentage = 20,
#   export = FALSE,
#   file_directory = NULL,
#   export_format = "csv"
# )

## ----eval=FALSE---------------------------------------------------------------
# search_pubmed(genes_list, terms_list, rank_method = "weighted")

## ----eval=FALSE---------------------------------------------------------------
# plot_heatmap(pubmed_search_results, file_directory = NULL, export = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# search_string_db(
#   genes_list,
#   species = 9606,
#   network_type = "full",
#   score_threshold = 400
# )

## ----eval=FALSE---------------------------------------------------------------
# plot_string_network(
#   string_db,
#   string_ids,
#   file_directory = NULL,
#   export = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# plot_clustering(string_results, file_directory = NULL, export = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# combine_summary(
#   pubmed_search_results,
#   string_results,
#   file_directory = NULL,
#   export_format = "csv",
#   export = FALSE,
#   threshold_percentage = 20
# )

## ----eval=FALSE---------------------------------------------------------------
# plot_connectivity_precedence(
#   combined_summary,
#   file_directory = NULL,
#   export = FALSE
# )

