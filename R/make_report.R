#' Automate knitr reporting
#' 
#' This function tries to automate some repetitive tasks in data analysis 
#' reporting. More precisely the function assumes that a differential expression
#' analysis has been performed for a series of contrasts, and that the results 
#' have been stored in a list, where each component's name is the contrast of 
#' interest. We may want to produce data tables and maybe some simple plots for
#' each contrasts, applying the same set of operations on each contrast. This 
#' function tries to achieve exactly this, generating an Rmd file that contains
#' headers reflecting the contrast under consideration, and producing data
#' tables for the results.
#' @param dataset List containing the outputs from \code{topTable} for several
#' contrasts of interest.
#' @param output Character string containing the full path to the file where
#' the Rmd report will be saved.
#' @param title Character string containing the title of the report.
#'
#' @return No return value, but an Rmd is saved in \code{output} as a 
#' side-effect.
#' @export
make_report <- function(dataset=NULL, output=NULL, title=NULL) {

  opener <- "```{r}"
  closer <- "```"
  src <- NULL
  
  # Header
  src <- paste("---",
               paste0("title: \"", title, "\""),      
               "author: Giovanni d'Ario",
               paste0("date: ", date()),
               "output:\n  html_document:\n    toc: true\n    toc_depth: 2",
               "---\n", sep = "\n")
    
#   src <- "# Summary of the Results\n"
  src <- paste(src, opener, 
               "library(knitr)", 
               "library(DT)",
               paste0("load(\"../", dataset, "\")"),
               closer, "\n", sep = "\n")

  dataset_name <- load(file = dataset)
  x <- get(dataset_name)
  nms <- names(x)
  
  for (i in seq_along(x)) {
    src <- paste(src, paste("## Contrast:", nms[i]), sep = "\n")
    src <- paste(
      src, opener, 
      paste0("datatable(", dataset_name, "[[", i, "]])"),
      closer, sep = "\n")
    src <- paste(src, "\n")
  }
  
  cat(src, file = output)
}
