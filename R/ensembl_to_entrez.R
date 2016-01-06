#' Reannotate a count matrix from Ensembl to Entrez ids
#' 
#' Take a count matrix generated in Ensembl space and 
#' re-annotate it in Entrez space.
#' 
#' This function takes a count matrix generated in Ensembl space
#' and returns an aggregated count matrix in Entrez space. The 
#' function relies on the mappings \code{human_mapping} and 
#' \code{mouse_mapping} to perform the aggregation. Note that this
#' implies a significant reduction in the number of features since,
#' in a typical experiment, count matrices can have > 60000 Ensembl
#' identifiers, while the mappings contain less than 30000.
#' Note also that this function aggregates using the mean. This 
#' happens because, in most cases, we need to go from Ensembl to 
#' Entrez \emph{after} an analysis is completed, i.e. when we are
#' dealing with a normalized and processed expression matrix.
#' When dealing with an actual count matrix, the right aggregation
#' would be a sum, not a mean.
#' @param x A count or an expression matrix produced in Ensembl 
#' space, with Ensembl gene ids on the rows.
#' @param species A character string specifying the species. As 
#' for now it can only be "human" or "mouse".
#' @param func A character string specifying the function to be 
#' applied to the matrix \code{x}. Possible values are \code{mean},
#' and \code{sum}. The latter is recommended for raw
#' count matrices from RNA-seq experiments.
#' @export
#' @author Giovanni d'Ario
#' @return A count matrix in Entrez space
ensembl_to_entrez <- function(x, species=c("human", "mouse"),
  func=c("sum", "mean")) {
  
  if (!is.matrix(x))
    stop("x must be a matrix with gene identifiers as rownames")
  
  species <- match.arg(species)
  func <- match.arg(func)
  
  mapping <- switch(species,
    human = human_mapping,
    mouse = mouse_mapping)
  
  dF <- data.frame(ensembl_gene_id = rownames(x), x,
    stringsAsFactors = FALSE)
  
  mapping_df <- dplyr::inner_join(mapping, dF)
  mapping_df$ensembl_gene_id <- NULL
  
  entrez_df <- dplyr::group_by(mapping_df, gene_id)
  
  entrez_df <- dplyr::summarise_each(entrez_df, 
    dplyr::funs_(func))
  
  out <- as.matrix(entrez_df[, -1])
  rownames(out) <- entrez_df$gene_id
  
  return(out)
}