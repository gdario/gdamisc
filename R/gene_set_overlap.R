##' Compute the proportion of overlap between two vectors
##'
##' This function takes two vectors and computes the proportion of
##' overlap as the ratio between the number of elements in the
##' intersection of the two objects, and the length of the smaller
##' vector.
##' @title Compute the proportion of overlapping elements in two
##' vectors.
##' @param x a vector of any class
##' @param y a vector of any class
##' @return a real number between 0 and 1 representing the fraction
##' of elements in the smaller vectors that are in common with the
##' larger vector.
##' @author Giovanni d'Ario
prop_overlap <- function(x, y) {
    l <- length(intersect(x, y))
    m <- min(length(x), length(y))
    return(l/m)
}

##' Compute the Jaccard index between two vectors
##'
##' Compute the Jaccard index between two vectors. The higher the
##' index, the more similar the two datasets.
##' @title Compute the Jaccard index
##' @param x vector of any type
##' @param y vector of any type
##' @return a real number between 0 and 1
##' @author Giovanni d'Ario
jaccard <- function(x, y) {
    i <- length(intersect(x, y))
    u <- length(union(x, y))
    return(i/u)
}

##' Compute the Sorensen Dice similarity coefficient
##'
##' @title Compute the Sorensen Dice similarity coefficient
##' @param x a vector of any type.
##' @param y a vector of any type.
##' @return a number between 0 and 1, where 1 means "more similar"
##' @author Giovanni d'Ario
dice <- function(x, y) {
    i <- length(intersect(x, y))
    m <- .5 * (length(x) + length(y))
    return(i/m)
}

##' Compute the overlap across gene sets
##'
##' This function takes a table of gene set IDs and their associated
##' gene IDs and compute the intersection across the different sets
##' @title Compute the overlap across two gene sets.
##' @param gsa_file the full path to the file containing the output of
##' the \code{gsa.py} analysis.
##' @param stat character string specifying the type of statistics to
##' be used for ranking the results of \code{gsa.py}.
##' @param similarity character, te type of similarity meausure to be
##' used to assess how similar two gene sets are in terms of content.
##' @param ranking character, the statistic used to rank the gene sets
##' in \code{gsa.py}. Default is the negative log_{10} of the p-value.
##' @param rank_by character, the type of rank used for ranking the
##' gene sets. This parameter only makes sense when using the (weighted)
##' Kolmogorov-Smirnov, which produces two statistics, Dm and Dp.
##' @param min_size the minimum size of a gene set.
##' @param max_size the maximum size of a gene set.
##' @param max_rank The largest rank of the output. This means that the
##' \code{igraph} object returned by the function will contain
##' information only about the gene sets up to the nth rank (not
##' necessarily n gene sets, in case of ties).
##' @param remove_zero_weight Logical. Should the edges with zero
##' weight be removed?
##' @return an object of class igraph.
##' @author Giovanni d'Ario
##' @export
gene_set_overlap <- function(gsa_file,
                             stat=c("KSw",
                                 "KS",
                                 "MannWhitneyU_U",
                                 "Wilcoxon_T",
                                 "tTest_t"),
                             similarity=c("dice",
                                 "jaccard",
                                 "overlap"),
                             ranking=c("signed_log10p",
                                 "fc"),
                             rank_by=c("rank_Dm",
                                 "rank_Dp",
                                 "rank"),
                             min_size=NULL,
                             max_size=NULL,
                             max_rank=20,
                             remove_zero_weight=TRUE) {

    ## Select the ranking criterion
    stat <- match.arg(stat)
    ranking <- match.arg(ranking)
    rank_by <- match.arg(rank_by)
    similarity <- match.arg(similarity)

    ## Read the gsa.py results and extract the columns of interest
    Data <- read_gsa(file = gsa_file,
                     stat = stat,
                     ranking = ranking,
                     min_size = min_size,
                     max_size = max_size)

    ## After filtering by min_size and max_size the ranking is no
    ## longer ordinal. We must fix this
    Data[[rank_by]] <- rank(Data[[rank_by]], ties.method = "min")

    ## Select only the highest ranking gene sets
    Data <- Data[order(Data[[rank_by]]), ]
    Data <- subset(Data, Data[[rank_by]] <= max_rank)

    ## Load the file containing the mapping between gene sets and
    ## Entrez Gene IDs
    gene_sets_file <- sub("\\.txt", "_set_genes.txt", gsa_file)
    Genes <- read.delim(file = gene_sets_file,
                        stringsAsFactors = FALSE)
    Genes <- subset(Genes, id_ %in% Data$id_)

    ## Compute the overlap/similarity across different gene sets
    gene_list <- split(Genes$gene_id, Genes$id_)

    ## Similarity function
    similarity_fun <- function(x, y, type) {
    	switch(type,
               jaccard = jaccard(x, y),
               dice = dice(x, y),
               overlap = prop_overlap(x, y))
    }

    ## Create a data frame containing the from and to gene set IDs
    ## and the value of the similarity score
    n <- length(gene_list)
    nms <- names(gene_list)

    similarity_df <- data.frame()

    for(i in 1:(n - 1)) {
    	for(j in (i + 1):n) {
            from <- nms[i]
            to <- nms[j]
            weight <- similarity_fun(x = gene_list[[i]],
                                     y = gene_list[[j]],
                                     type = similarity)

            similarity_df <- rbind(similarity_df,
                                   data.frame(from, to, weight,
                                              stringsAsFactors = FALSE))
    	}
    }

    ## Store the gene set ID and the gene set name
    gene_set_info <- unique(subset(Data,
                                   select = c("id_",
                                       rank_by,
                                       "name",
                                       "gene_set_size")))

    gene_set_info$id_ <- as.character(gene_set_info$id_)

    g <- graph.edgelist(el = as.matrix(
                            similarity_df[c('from', 'to')]),
                        directed = FALSE)

    E(g)$weight <- similarity_df$weight

    idx <- match(V(g)$name, gene_set_info$id_)

    if(any(is.na(idx)))
    	stop("There is a problem with the identifiers")

    gene_set_info <- gene_set_info[idx, ]

    V(g)$description <- gene_set_info$name
    V(g)$gene_set_size <- gene_set_info$gene_set_size
    V(g)$rank <- gene_set_info[[rank_by]]

    ## Remove the edges with a zero weight
    if(remove_zero_weight) {
        idx <- which(E(g)$weight == 0)
        g <- g - edges(idx)
    }

    return(g)
}
