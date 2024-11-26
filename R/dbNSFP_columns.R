#' List Available Columns in dbNSFP
#'
#' This function lists all available column names in a dbNSFP file and standardizes them.
#'
#' @param dbnsfp_file Path to the bgzipped dbNSFP file.
#'
#' @return A character vector of column names, with standardization for specific cases (e.g., `#chr` to `chr`).
#' @export
dbNSFP_columns <- function(dbnsfp_file) {
  
  # Open the Tabix file
  tabix_file <- TabixFile(dbnsfp_file)
  
  # Read the header
  header <- headerTabix(tabix_file)
  
  # Extract the column names
  column_names <- strsplit(header$header, "\t")[[1]]
  
  # Standardize the column names (e.g., remove `#` from `#chr`)
  column_names <- gsub("^#chr$", "chr", column_names)
  
  return(column_names)
}
