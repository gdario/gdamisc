#' Ensembl-Entrez Id mapping for mouse
#' 
#' A data frame containing the mapping between Ensembl gene
#' ids and Entrez ids for mus musculus.
#' 
#' This data frame contains a mapping between the Ensembl gene 
#' ids and the Entrez gene ids. It is used by the 
#' \code{ensembl_to_entrez} function, that takes a count matrix
#' based on Ensembl and returns it aggregating the Entrez ids
#' when necessary. The data have been extracted from Hithub.
#' 
#' @format a data frame with 20313 rows and 2 columns:
#' \describe{
#'  \item{ensembl_gene_id}{The Ensemble gene identifier.}
#'  \item{gene_id}{The Entrez gene identifier}
#' }
"mouse_mapping"