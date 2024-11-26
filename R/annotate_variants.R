#' Annotate Variants Using dbNSFP
#'
#' This function annotates a set of variants using the dbNSFP database. 
#' It supports queries either in HGVSg format or as separate columns.
#'
#' @param query A data frame or vector. If `is_HGVSg = TRUE`, provide a vector of HGVSg strings. 
#'   Otherwise, provide a data frame with columns: `chr`, `start`, `end`, `ref`, `alt`.
#' @param dbnsfp_file Path to the bgzipped and indexed dbNSFP file.
#' @param columns A character vector of column names from dbNSFP to include in the annotation. Default is NULL, which uses all available columns.
#' @param is_HGVSg Boolean. If `TRUE`, treats input as HGVSg strings. Default is `FALSE`.
#' @param chunk_size Number of variants to process per chunk. Default is 1000.
#' @param workers Number of cores to use for parallelization. Default is 6.
#'
#' @return A data frame with the original query data and the selected dbNSFP annotations.
#' @export
annotate_variants <- function(query, dbnsfp_file, columns = NULL, is_HGVSg = FALSE, chunk_size = 1000, workers = 6) {
  
  # Parse HGVSg if specified
  if (is_HGVSg) {
    message("HGVSg mode detected. Parsing HGVSg strings...")
    query_data <- parse_HGVSg(query)
  } else {
    query_data <- query
    if (!all(c("chr", "start", "end", "ref", "alt") %in% colnames(query_data))) {
      stop("Input data must contain columns: chr, start, end, ref, alt when is_HGVSg = FALSE.")
    }
  }
  
  # Check if the dbNSFP file exists and is valid
  if (!file.exists(dbnsfp_file)) {
    stop("The specified dbNSFP file does not exist. Please provide a valid file path.")
  }
  index_file <- paste0(dbnsfp_file, ".tbi")
  if (!file.exists(index_file)) {
    stop("The specified dbNSFP file does not have a valid Tabix index (.tbi file). Please index the file using tabix.")
  }
  
  # Get available columns from dbNSFP
  available_columns <- dbNSFP_columns(dbnsfp_file)
  
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
  
  # Split data into chunks (adjusted for vectors)
  if (is.vector(query_data)) {
    chunk_indices <- split(seq_along(query_data), (seq_along(query_data) - 1) %/% chunk_size)
  } else {
    chunk_indices <- split(seq_len(nrow(query_data)), (seq_len(nrow(query_data)) - 1) %/% chunk_size)
  }
  
  # Define query function
  query_dbnsfp <- function(indices) {
    result_list <- vector("list", length(indices))
    for (i in seq_along(indices)) {
      tryCatch({
        index <- indices[i]
        region <- GRanges(
          seqnames = query_data$chr[index],
          ranges = IRanges(
            start = query_data$start[index],
            end = query_data$end[index]
          )
        )
        result <- scanTabix(TabixFile(dbnsfp_file), param = region)
        
        if (length(result) > 0) {
          # Parse result and select specified columns
          parsed_result <- strsplit(unlist(result), "\t")[[1]]
          names(parsed_result) <- available_columns
          
          # Ensure only the requested columns are selected
          selected_result <- tryCatch(
            parsed_result[columns],
            error = function(e) {
              # Handle subscript out of bounds error
              message(sprintf("Subscript out of bounds for variant %d: %s", index, e$message))
              setNames(rep(NA, length(columns)), columns)
            }
          )
          result_list[[i]] <- selected_result
        } else {
          # If no matching entry is found, return NA for all selected columns
          message(sprintf("No match found for variant %d in dbNSFP.", index))
          result_list[[i]] <- setNames(rep(NA, length(columns)), columns)
        }
      }, error = function(e) {
        # General error handling
        message(sprintf("Error querying variant %d: %s", indices[i], e$message))
        result_list[[i]] <- setNames(rep(NA, length(columns)), columns)
      })
    }
    return(do.call(rbind, result_list))
  }
  
  # Process chunks in parallel
  all_annotations <- bplapply(chunk_indices, query_dbnsfp)
  
  # Combine all chunks into a single data frame
  annotation_data <- do.call(rbind, all_annotations)
  
  # Merge annotations with original data
  result <- cbind(query_data, annotation_data)
  return(result)
}
