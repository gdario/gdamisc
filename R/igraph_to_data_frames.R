#' This function onverts an igraph objects into a pair of data
#' frames, one for the edges and one for the vertices. The user
#' provides an output file name and the function creates two files
#' with the same basename but appending _edges and _vertices it.
#'
#' Put the details here
#' @title Convert an igraph objects into a pair of data frames
#' @param graph an object of class \code{igraph}
#' @param output_file
igraph_to_data_frames <- function(graph, output_file) {
	el <- get.edgelist(graph, names = TRUE)
	va <- vertex.attributes(graph)
	ea <- edge.attributes(graph)
	edf <- data.frame(from = el[, 1], to = el[, 2], ea)
	edge_file <- paste0(output_file, "_edges.txt")
	vertex_file <- paste0(output_file, "_vertices.txt")
	write.table(edf, file = edge_file, sep = "\t", quote = FALSE,
				row.names = FALSE)
	write.table(va, file = vertex_file, sep = "\t", quote = FALSE,
				row.names = FALSE)
}
