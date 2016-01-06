#' Convert a list into a GMT file
#'
#' This function converts a list into a GMT file
#' @param x the list to be converted into a GMT file.
#' @param gmt_file the file where the GMT file should be saved
#' @title Convert a list into a GMT file.
#' 
#' **TO DO** Add an option to select if one or two elements
#' should be considered as identifiers at the beginning of each
#' list component.
#' @seealso \code{gmt_to_list}
#' @author Giovanni d'Ario
#' @export
#' @examples
#' # gmt_file = "" prints on the standard output.
#' gene_list <- list(c("process_A", "gene_1", "gene_2"),
#'                   c("process_B", "gene_3", "gene_4"))
#' list_to_gmt(gene_list, gmt_file = "")
list_to_gmt <- function(x, gmt_file) {

	cat(file = gmt_file, paste(x[[1]], sep = "",
							   collapse = "\t"),
		"\n")
	for(i in 2:length(x)) {
		cat(file = gmt_file, paste(x[[i]], sep = "",
								   collapse = "\t"),
			"\n", append = TRUE)
	}
}
