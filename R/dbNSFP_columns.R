#' List Available Columns in dbNSFP
#'
#' This function lists all available column names in a dbNSFP file.
#'
#' @param dbnsfp_file Path to the bgzipped dbNSFP file.
#'
#' @return A character vector of column names if the file is valid and has a proper header.
#' @export
dbNSFP_columns <- function(dbnsfp_file) {
  
  # Check if file exists
  if (!file.exists(dbnsfp_file)) {
    stop("The specified dbNSFP file does not exist. Please provide a valid file path.")
  }
  
  # Check if the file has a valid Tabix index
  index_file <- paste0(dbnsfp_file, ".tbi")
  if (!file.exists(index_file)) {
    stop("The specified dbNSFP file does not have a valid Tabix index (.tbi file). Please index the file using tabix.")
  }
  
  # Try opening the Tabix file
  tryCatch({
    tabix_file <- TabixFile(dbnsfp_file)
  }, error = function(e) {
    stop("Failed to open the dbNSFP file. Ensure it is a valid bgzipped file.")
  })
  
  # Try reading the header
  header <- tryCatch({
    headerTabix(tabix_file)
  }, error = function(e) {
    stop("Failed to read the header of the dbNSFP file. Ensure the file is correctly formatted.")
  })
  
  # Check if the header is non-empty
  if (length(header) == 0) {
    stop("The dbNSFP file appears to have an empty or missing header. Please check the file format.")
  }
  
  # Extract the column names
  column_names <- tryCatch({
    strsplit(header$header, "\t")[[1]]
  }, error = function(e) {
    stop("Failed to parse the header of the dbNSFP file. Ensure the header line is tab-delimited.")
  })
  
  # Check if any columns were extracted
  if (length(column_names) == 0) {
    stop("No column names could be extracted from the dbNSFP file. Please check the file format.")
  }
  
  # Return the column names
  return(column_names)
}
