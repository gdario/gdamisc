test_that("Gene identifiers are correctly mapped via biomaRt", {
  set.seed(123)
  
  ensg <- c("ENSG00000197421", "ENSG00000140323", "ENSG00000152213", 
           "ENSG00000278558", "ENSG00000125775")
  
  # Create a matrix version of the data
  m <- matrix(rnorm(25), nrow = 5)
  colnames(m) <- paste0("S", 1:5)
  
  # Create a data frame version of the data
  d <- data.frame(ensembl_gene_id = ensg, m)
  rownames(m) <- ensg
  
  # Generate the output for different combinations of input class and options.
  out_m <- map_identifiers(x = m, collapse = FALSE)
  out_d <- map_identifiers(x = d, collapse = FALSE)
  out_m_c <- map_identifiers(x = m, collapse = TRUE)
  out_d_c <- map_identifiers(x = d, collapse = TRUE)
  
  expect_that(out_m, is_a("matrix"))
  expect_that(out_d, is_a("matrix"))
  expect_that(out_m_c, is_a("matrix"))
  expect_that(out_d_c, is_a("matrix"))
  
  expect_that(dim(out_m), equals(c(5, 5)))
  expect_that(dim(out_d), equals(c(5, 5)))
  expect_that(dim(out_m_c), equals(c(5, 5)))
  expect_that(dim(out_d_c), equals(c(5, 5)))
})