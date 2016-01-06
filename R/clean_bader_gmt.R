#' Convert the GMT from the Bader's lab into a format that can be
#' read and processed by \code{gsa.py}.
#'
#' The Bader's lab (\url{http://baderlab.org/GeneSets}) contains a
#' large set of GMT files that have been mapped to mouse. However,
#' when run in \code{gsa.py}, these files cause the program to 
#' crash, most probably due to the presence of non UTF-8 encodings.
#' This program reads a Bader's GMT file, cleans it removing some
#' of the redundant (and lengthy) strings, checking the encoding,
#' and saves the cleaned version.
#' @param input_file path to the "dirty" GMT file.
#' @param output_file where the "cleaned" output file should be
#' saved.
#' @title Clean Bader's GMT files
#' @author Giovanni d'Ario
#' @export
clean_bader_gmt <- function(input_file, output_file) {
	
	gmt <- scan(input_file, what = "", sep = "\n")
        gmt <- lapply(gmt, iconv, from = "LATIN1", to = "UTF-8")
	gmt <- sapply(gmt, strsplit, split = "\t")
	
	out <- lapply(gmt, function(x) {
		y <- x
		y[1] <- toupper(y[2])
		return(y)
	})

	## Save the cleaned gmt file
	list_to_gmt(x = out, gmt_file = output_file)
}
