test_that("annotate_variants works as expected", {
  query_data <- data.frame(
    chr = c("1", "2"),
    start = c(12345, 67890),
    end = c(12345, 67890),
    ref = c("A", "G"),
    alt = c("T", "C")
  )

  # Replace with a small test dbNSFP file
  dbnsfp_file <- "path/to/test/dbNSFP.gz"

  result <- annotate_variants(query_data, dbnsfp_file, chunk_size = 1, workers = 2)

  expect_true(nrow(result) == nrow(query_data))
})
