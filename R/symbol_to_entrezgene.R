##' Convert a GMT from symbol-based to Entrez-based
##'
##' This function reads a GMT file and uses biomaRt to convert the
##' gene symbols (MGI for mouse and HGNC for human) to the
##' corresponding Entrez IDs. It returns a list that can be then
##' saved using \code{list_to_gmt} into a new GMT file
##' @title Convert a GMT from symbol-based to Entrez-based
##' @param gmt_file full path to a GMT file
##' @param species character, either "mouse" (default) or "human"
##' @return a list containing the Entrez IDs of the GMT file
##' provided in input
##' @author Giovanni d'Ario
##' @export
symbol_to_entrezgene <- function(gmt_file,
                                 species=c("mouse", "human")) {
    require(biomaRt)

    species <- match.arg(species)

    message("connecting to biomart")
    if(species == "mouse") {
        dataset <- "mmusculus_gene_ensembl"
        symbol <- "mgi_symbol"
    }
    if(species == "human") {
        dataset <- "hsapiens_gene_ensembl"
        symbol <- "hgnc_symbol"
    }
    
    ensembl <- useMart(biomart = "ensembl",
                       dataset = dataset)
    gmt_list <- gmt_to_list(gmt_file)

    ## I have found gmt files having gene lists with no genes!
    ll <- sapply(gmt_list, length)
    gmt_list <- gmt_list[ll > 2]
    
    ## Turn the gene symbols into
    message("Converting gene symbols to Entrez ids")
    out_list <- lapply(gmt_list, function(x) {
        y <- x[-c(1,2)]
        eid <- getBM(attributes = "entrezgene",
                     filters = symbol,
                     values = y,
                     mart = ensembl)
        ## Repeat the second element of the GMT in the first
        ## position to avoid problems when reading it.
        out <- c(x[2], x[2], unique(eid$entrezgene))
        return(out)
    })
    return(out_list)
}
