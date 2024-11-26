#' Parse HGVSg Strings
#'
#' This function parses HGVSg-formatted strings into separate components: 
#' chromosome, start, end, reference allele, and alternate allele.
#'
#' @param hgvs_strings A character vector of HGVSg strings (e.g., `chr1:g.12345A>T`, `chr12:g.55723908del`).
#'
#' @return A data frame with columns: `chr`, `start`, `end`, `ref`, `alt`. Invalid entries will have NA values.
#' @export
parse_HGVSg <- function(hgvs_strings) {
  # Regular expressions for different types of HGVSg strings
  substitution_pattern <- "^chr(\\w+):g\\.(\\d+)([ACGTN])>([ACGTN])$"
  deletion_pattern <- "^chr(\\w+):g\\.(\\d+)(_\\d+)?del$"
  duplication_pattern <- "^chr(\\w+):g\\.(\\d+)(_\\d+)?dup$"
  insertion_pattern <- "^chr(\\w+):g\\.(\\d+)_\\d+ins([ACGTN]+)$"
  
  parsed_data <- lapply(hgvs_strings, function(hgvs) {
    if (grepl(substitution_pattern, hgvs)) {
      # Substitution
      matches <- regmatches(hgvs, regexec(substitution_pattern, hgvs))[[1]]
      return(data.frame(chr = matches[2], start = as.numeric(matches[3]), 
                        end = as.numeric(matches[3]), ref = matches[4], alt = matches[5]))
    } else if (grepl(deletion_pattern, hgvs)) {
      # Deletion
      matches <- regmatches(hgvs, regexec(deletion_pattern, hgvs))[[1]]
      start <- as.numeric(matches[3])
      end <- ifelse(matches[4] != "", as.numeric(sub("_(\\d+)", "\\1", matches[4])), start)
      return(data.frame(chr = matches[2], start = start, end = end, ref = "del", alt = ""))
    } else if (grepl(duplication_pattern, hgvs)) {
      # Duplication
      matches <- regmatches(hgvs, regexec(duplication_pattern, hgvs))[[1]]
      start <- as.numeric(matches[3])
      end <- ifelse(matches[4] != "", as.numeric(sub("_(\\d+)", "\\1", matches[4])), start)
      return(data.frame(chr = matches[2], start = start, end = end, ref = "", alt = "dup"))
    } else if (grepl(insertion_pattern, hgvs)) {
      # Insertion
      matches <- regmatches(hgvs, regexec(insertion_pattern, hgvs))[[1]]
      return(data.frame(chr = matches[2], start = as.numeric(matches[3]), 
                        end = as.numeric(matches[3]), ref = "", alt = matches[4]))
    } else {
      # Invalid format
      return(data.frame(chr = NA, start = NA, end = NA, ref = NA, alt = NA))
    }
  })
  
  # Combine parsed results into a data frame
  parsed_df <- do.call(rbind, parsed_data)
  rownames(parsed_df) <- NULL
  invalid_rows <- which(is.na(parsed_df$chr))
  if (length(invalid_rows) > 0) {
    warning(sprintf(
      "Invalid HGVSg format detected for entries: %s. Ensure all strings follow the expected formats.",
      paste(hgvs_strings[invalid_rows], collapse = ", ")
    ))
  }
  return(parsed_df)
}
