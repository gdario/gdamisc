##' Compute the ranking of the gene sets produced by gsa.py
##'
##' This function computes the ranks of the gene sets loaded by
##' the readGsaPy function. The ranking is simply produced by 
##' separately ranking the statistics and its associated 
##' -log10(adjusted p-value) and by then
##' ranking the sum of the two ranks.
##' @title Compute the ranking of the gene sets produced by gsa.py
##' @param x a dataset loaded by readGsaPy
##' @param stat the statistics of interest.
##' @return an integer vector indicating the ranks (1 = highest ranking)
##' of the gene sets.
##' @author Giovanni d'Ario
##' @export
compute_rank <- function(x, stat) {

     idxPv <- grep(paste0(stat, "_pval_adj_neglog10$"),
                  colnames(x))

#     idxPv <- grep(paste0(stat, "_pval_neglog10$"),
#                   colnames(x))

    ## We want to rank in decreasing order, but this is not
    ## available in R, therefore we rank the negative of the
    ## values.
    rankPv <- rank(-x[[idxPv]], ties.method = "min")

    idxStat <- grep(paste0(stat, "$"), colnames(x))
    rankStat <- rank(-x[[idxStat]], ties.method = "min")

    finalRank <- rank(rankPv + rankStat, ties.method = "min")

    return(finalRank)
}

##' Add the ranking to the GSA.py data
##'
##' This function
##' @title Add the ranking columns to the gsa.py output
##' @param x a data frame containing the output of gsa.py read
##' by \code{read_gsa}
##' @param stat character string. The statistics used for the
##' ranking. Defaults to "KSw".
##' @return a data frame identical to the one provided in input
##' but with one or two additional columns containing the ranking.
##' @author Giovanni d'Ario
add_ranking <- function(x, stat) {
    ## If the statistic is KSw or KS you must treat the Dp and the Dm
    ## cases separately.

    if((stat == "KS") | (stat == "KSw")) {
        rankDm <- compute_rank(x, stat = paste0(stat, "_Dm"))
        rankDp <- compute_rank(x, stat = paste0(stat, "_Dp"))
        out <- cbind(rank_Dp = rankDp, rank_Dm = rankDm, x)
    } else {
        rank <- compute_rank(x, stat = stat)
        out <- cbind(rank, x)
    }

    return(out)
}

##' Read into R the output of gsa.py
##'
##' This function takes the full path to a gsa.py output file and
##' returns a data frame containing the results of the analysis,
##' selecting a subset of the columns. Please note
##' that this function is limited to two-group comparisons. Extensions
##' to more complex designs will be added in the future.
##' @title Read the GSA.py output files into R and add the ranking
##' @param file full path to the gsa.py output file.
##' @param stat character string, the statistical test used by the
##' gsa.py script that should be used for the ranking. The default is
##' the weighted Kolmogorov Smirnov.
##' @param ranking the criterion that should be used to rank gene sets.
##' The default is "signed_log10p".
##' @param min_size optional. The minimum size of the gene sets that should
##' be included in the output.
##' @param max_size optional. The maximum size of the gene sets that should
##' be included in the output.
##' @return a data frame restricted to a subset of the columns present
##' in the gsa.py output file,
##' @author Giovanni d'Ario
read_gsa <- function(file,
					 stat=c("KSw",
					 	   "KS",
					 	   "MannWhitneyU_U",
					 	   "Wilcoxon_T",
					 	   "tTest_t"),
					 ranking=c("fc", "signed_log10p"),
					 min_size=NULL,
					 max_size=NULL) {

    stat <- match.arg(stat)
    ranking <- match.arg(ranking)
    ranking <- paste0("^", ranking)

    nms <- read.delim(file, header = TRUE, nrows = 1)
    nms <- names(nms)
    Data <- read.delim(file, header = FALSE, skip = 2,
                       stringsAsFactors = FALSE)
    names(Data) <- nms

    cols <- c("id_", "category",
              "name", "xref",
              "gene_set_size")

    idxCol1 <- match(cols, colnames(Data))
    idxCol2 <- grep(paste(ranking, stat, sep = "_"),
                    names(Data))

    out <- Data[c(idxCol1, idxCol2)]

    ## Remove the columns with the unadjusted p-value
    idx <- grep("pval_neglog10$", names(out))
    out <- out[, -idx]

    ## Add the ranking
    out <- add_ranking(out, stat = stat)

    ## Filter by gene set size
    if(!is.null(min_size))
        out <- subset(out, gene_set_size >= min_size)
    if(!is.null(max_size))
        out <- subset(out, gene_set_size <= max_size)

    return(out)
}
