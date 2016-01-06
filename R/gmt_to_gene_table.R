##' Transform a GMT file into a rectangular matrix
##'
##' This function can receive either a character string
##' representing the full path to a GMT file (that will be then
##' read by \code{gmt_to_list}, or a list (obtained by running
##' \code{gmt_to_list}. It then creates a rectangular matrix
##' where the fist \code{n_descriptors} columns contain the
##' descriptors of the gene set (for example the name and the
##' description), and the remaining column contains the gene ID
##' of the gene in the gene set. 
##' @title Transform a GMT file into a rectangular matrix
##' @param x either the path to a GMT file, or a list generated
##' by a call to gmt_to_list.
##' @param n_descriptors integer, the number of descriptor
##' column present in the GMT file.
##' @param names A character vector containing the names that
##' should appear in the final table. If NULL, these will be
##' automatically generated as "descriptor_1", "descriptor_2",
##' ..., "gene_id"
##' @return a character matrix containing the descriptors and
##' the gene ids for each gene in each gene set.
##' @author Giovanni d'Ario
##' @export
gmt_to_gene_table <- function(x=NULL,
                              n_descriptors=NULL,
                              names=NULL,
                              gene_id_first=TRUE) {

    ## The argument x can either be the path to a gmt file
    ## or a list obtained by running gmt_to_list
    if(is.character(x))
        gmt_list <- gmt_to_list(gmt_file = x)
    if(is.list(x))
        gmt_list <- x
    
    if(is.null(n_descriptors))
        stop("n_descriptors must indicate the number of non-gene id columns")
    
    f <- function(x) {
        ll <- length(x) - n_descriptors
        descriptors <- x[1:n_descriptors]
        gene_ids <- x[-c(1:n_descriptors)]
        dd <- rep(descriptors, each = ll)
        dd <- matrix(dd, ncol = n_descriptors)
        out <- cbind(dd, gene_ids)
        out
    }

    gene_table <- lapply(gmt_list, f)
    res <- do.call("rbind", gene_table)

    if(is.null(names))
        colnames(res) <- c(paste("descriptor",
                              1:n_descriptors, sep = "_"),
                        "gene_id")

    if (gene_id_first) {
        res <- cbind(res[, n_descriptors + 1],
                     res[, 1:n_descriptors])

        colnames(res) <- c("gene_id",
                           paste("descriptor",
                              1:n_descriptors, sep = "_"))
    }
    res
}
