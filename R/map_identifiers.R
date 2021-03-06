#' Map between two gene identifiers using Biomart
#'
#' This function takes an expression matrix or a data frame containing some
#' numerical values and some identifiers. If the input \code{x} is a matrix the
#' identifiers are expected to be the row names of the matrix. If \code{x} is 
#' a data frame, the function expects that \code{x} contains a column whose 
#' name is the same as the string contained in \code{id1}, the name of the 
#' identifiers we want to move away from (default is \code{ensembl_gene_id}).
#' Optionally, the function collapses the numerical values applying the 
#' function specified in \code{func} to the identifiers of type \code{id2} that
#' appear multiple times. For RNA-Seq raw data the best function to use is
#' \code{sum}. For normalized and log-scale microarray data one can use either
#' \code{mean} or \code{median}.
#' 
#' @param x a matrix or a data frame containing the expression values.
#' @param id1 Character string containing the name of the first gene
#' identifier. Defaults to \code{ensembl_gene_id}.
#' @param id2 Character string containing the name of the second gene
#' identifier. Defaults to \code{hgnc_symbol}.
#' @param values character string containing the \code{id1} identifiers. If
#' not specified, and if \code{x} is a matrix, they will be assumed to be the
#' row names of the matrix. If \code{x} is a data frame, the function assumes 
#' that it contains a column whose name is the same as \code{id1}.
#' @param version Integer specifying the Ensembl version number. Defaults to 83.
#' @param organism Character string specifying the organism. Possibile values
#' are \code{hsapiens} and \code{mmusculus}. Defaults to \code{hsapiens}.
#' @param collapse Logical. Should the rows of x be collapsed with respect
#' to \code{id2}?
#' @param func Function. The function used to collapse the rows of x. 
#' Defaults to \code{sum}.
#' @param as_matrix Logical. Should the output be a matrix?
#' @param rm_blank_ids Logical. Should the blank gene IDs be removed? Defaults
#' to \code{TRUE}.
#' @return a data frame containing the (possibly collapsed) expression matrix
#' mapped to the \code{id2} identifiers.
#' @export
#' @examples
#' ensg <- c("ENSG00000197421", "ENSG00000140323", "ENSG00000152213", 
#' "ENSG00000278558", "ENSG00000125775")
#' # Create an expression matrix...
#' m <- matrix(rnorm(25), nrow = 5)
#' colnames(m) <- paste0("S", 1:5)
#' # ...and a data frame.
#' d <- data.frame(ensembl_gene_id = ensg, m)
#' rownames(m) <- ensg
#' out <- map_identifiers(x = m, collapse = FALSE)

map_identifiers <- function(x=NULL,
                            id1='ensembl_gene_id', 
                            id2='hgnc_symbol', 
                            values=NULL,
                            version=83, 
                            organism=c('hsapiens', 'mmusculus'),
                            collapse=TRUE,
                            func=sum,
                            as_matrix=TRUE,
                            rm_blank_ids=TRUE) {

  require(magrittr)
  organism <- match.arg(organism)
  
  if (is.null(values)) {
    if (is.matrix(x)) {
      values <- rownames(x)
    } else if (is.data.frame(x)) {
      values <- as.character(x[[id1]])
    } else {
      stop("Which identifiers shall I map?")
    }
  }
  
  # Load the organism-specific biomart
  ensembl <- biomaRt::useEnsembl(
    biomart = 'ensembl', 
    dataset = paste0(organism, '_gene_ensembl'),
    version = version
  )
  
  # Map the identifiers using biomaRt
  id1_to_id2 <- biomaRt::getBM(attributes = c(id1, id2), 
                               filters = id1,
                               values = values, 
                               mart = ensembl)
  
  if (is.matrix(x)) {
    dF <- data.frame(tmp = rownames(x), x)
    names(dF) <- c(id1, colnames(x))
  } else {
    dF <- x
  }
  
  dF <- unique(dplyr::inner_join(id1_to_id2, dF))
  
  # Remove the id1 identifiers
  dF <- dplyr::select(dF, -matches(id1))
  
  # If required, collapse the identifiers having the same symbol by using 
  # the function speicified in `func`. Sums are suitable for RNA-Seq data.
  if (collapse) {
    dF <- dF %>% 
      dplyr::group_by_(id2) %>% 
      dplyr::summarise_each(dplyr::funs(func)) %>%
      dplyr::ungroup()
  }
  
  if (rm_blank_ids) {
    idx <- dF[[id2]] != ""
    dF <- dF[idx, ]
  }
  
  if (as_matrix) {
    rownames(dF) <- dF[[id2]]
    dF <- dplyr::select(dF, -matches(id2))
    dF <- as.matrix(dF)
  }
  
  return(dF)
}
