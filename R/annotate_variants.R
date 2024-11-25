#' Annotate Variants Using dbNSFP
#'
#' This function annotates a set of variants using the dbNSFP database.
#'
#' @param query_data A data frame containing variant information with columns:
#'   `chr` (chromosome), `start` (start position), `end` (end position), `ref` (reference allele), and `alt` (alternate allele).
#' @param dbnsfp_file Path to the bgzipped and indexed dbNSFP file.
#' @param chunk_size Number of variants to process per chunk. Default is 1000.
#' @param workers Number of cores to use for parallelization. Default is 6.
#'
#' @return A data frame with the original query data and dbNSFP annotations.
#' @export
annotate_variants <- function(query_data, dbnsfp_file, chunk_size = 1000, workers = 6) {
  library(VariantAnnotation)
  library(Rsamtools)
  library(BiocParallel)

  # Register parallel backend
  param <- MulticoreParam(workers = workers)
  register(param)

  # Split data into chunks
  chunks <- split(query_data, (seq_len(nrow(query_data)) - 1) %/% chunk_size)

  # Ensure Tabix file is accessible
  tabix_file <- TabixFile(dbnsfp_file)

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
        result <- scanTabix(tabix_file, param = region)
        if (length(result) > 0) {
          result_list[[i]] <- unlist(result)
        } else {
          result_list[[i]] <- NA
        }
      }, error = function(e) {
        message(sprintf("Error querying variant %d: %s", i, e$message))
        result_list[[i]] <- NA
      })
    }
    return(result_list)
  }

  # Process chunks in parallel
  process_chunk <- function(chunk) {
    annotations <- query_dbnsfp(chunk)
    annotation_data <- do.call(rbind, lapply(annotations, function(x) {
      if (!is.na(x)) {
        strsplit(x, "\t")[[1]]
      } else {
        rep(NA, 1)  # Adjust length to match dbNSFP fields
      }
    }))
    return(annotation_data)
  }

  # Run annotation
  all_annotations <- bplapply(chunks, process_chunk)
  annotation_data <- do.call(rbind, all_annotations)

  # Merge annotations with original data
  result <- cbind(query_data, annotation_data)
  return(result)
}
