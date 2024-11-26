#' Annotate Variants Using dbNSFP
#'
#' This function annotates a set of variants using the dbNSFP database. 
#' It allows the user to specify which columns from the dbNSFP database 
#' should be used for annotation.
#'
#' @param query A data frame containing variant information with columns: 
#'   `chr` (chromosome), `start` (start position), `end` (end position), `ref` (reference allele), and `alt` (alternate allele).
#' @param dbnsfp_file Path to the bgzipped and indexed dbNSFP file.
#' @param columns A character vector of column names from dbNSFP to include in the annotation. Default is NULL, which uses all available columns.
#' @param chunk_size Number of variants to process per chunk. Default is 1000.
#' @param workers Number of cores to use for parallelization. Default is 6.
#'
#' @return A data frame with the original query data and the selected dbNSFP annotations.
#' @export
annotate_variants <- function(query, dbnsfp_file, columns = NULL, chunk_size = 1000, workers = 6) {
  library(VariantAnnotation)
  library(Rsamtools)
  library(BiocParallel)
  
  # Check if the dbNSFP file exists and is valid
  if (!file.exists(dbnsfp_file)) {
    stop("The specified dbNSFP file does not exist. Please provide a valid file path.")
  }
  index_file <- paste0(dbnsfp_file, ".tbi")
  if (!file.exists(index_file)) {
    stop("The specified dbNSFP file does not have a valid Tabix index (.tbi file). Please index the file using tabix.")
  }
  
  # Get available columns from dbNSFP
  available_columns <- list_dbNSFP_columns(dbnsfp_file)
  
  # Validate selected columns
  if (!is.null(columns)) {
    missing_columns <- setdiff(columns, available_columns)
    if (length(missing_columns) > 0) {
      stop(sprintf(
        "The following columns are missing from the dbNSFP file: %s", 
        paste(missing_columns, collapse = ", ")
      ))
    }
  } else {
    # If no columns are specified, use all available columns
    columns <- available_columns
  }
  
  # Register parallel backend
  param <- MulticoreParam(workers = workers)
  register(param)
  
  # Split data into chunks
  chunks <- split(query, (seq_len(nrow(query)) - 1) %/% chunk_size)
  
  # Define query function
  query_dbnsfp <- function(query_chunk) {
    result_list <- vector("list", nrow(query_chunk))
    for (i in seq_len(nrow(query_chunk))) {
      tryCatch({
        region <- GRanges(
          seqnames = query_chunk$chr[i],
          ranges = IRanges(
            start = query_chunk$start[i],
            end = query_chunk$end[i]
          )
        )
        result <- scanTabix(TabixFile(dbnsfp_file), param = region)
        
        if (length(result) > 0) {
          # Parse result and select specified columns
          parsed_result <- strsplit(unlist(result), "\t")[[1]]
          names(parsed_result) <- available_columns
          selected_result <- parsed_result[columns]
          result_list[[i]] <- selected_result
        } else {
          result_list[[i]] <- setNames(rep(NA, length(columns)), columns)
        }
      }, error = function(e) {
        message(sprintf("Error querying variant %d: %s", i, e$message))
        result_list[[i]] <- setNames(rep(NA, length(columns)), columns)
      })
    }
    return(do.call(rbind, result_list))
  }
  
  # Process chunks in parallel
  all_annotations <- bplapply(chunks, query_dbnsfp)
  
  # Combine all chunks into a single data frame
  annotation_data <- do.call(rbind, all_annotations)
  
  # Merge annotations with original data
  result <- cbind(query, annotation_data)
  return(result)
}
