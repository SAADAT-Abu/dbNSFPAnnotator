#' Parse HGVSg Strings
#'
#' This function parses HGVSg-formatted strings into separate components: 
#' chromosome, start, end, reference allele, and alternate allele.
#'
#' @param query A character vector of HGVSg strings (e.g., `chr1:g.12345A>T`).
#'
#' @return A data frame with columns: `chr`, `start`, `end`, `ref`, `alt`.
#' @export
parse_HGVSg <- function(query) {
  # Regular expression to extract components
  hgvs_pattern <- "^chr(\\w+):g\\.(\\d+)([ACGTN])>([ACGTN])$"
  
  # Parse HGVSg strings
  parsed <- regmatches(query, regexec(hgvs_pattern, query))
  
  # Check for invalid HGVSg strings
  if (any(sapply(parsed, length) != 5)) {
    invalid_entries <- which(sapply(parsed, length) != 5)
    stop(sprintf(
      "Invalid HGVSg format detected for entries: %s. Ensure all strings follow the format: 'chr1:g.12345A>T'.",
      paste(query[invalid_entries], collapse = ", ")
    ))
  }
  
  # Create data frame from parsed components
  data.frame(
    chr = sapply(parsed, `[[`, 2),
    start = as.numeric(sapply(parsed, `[[`, 3)),
    end = as.numeric(sapply(parsed, `[[`, 3)),  # Single position variants
    ref = sapply(parsed, `[[`, 4),
    alt = sapply(parsed, `[[`, 5)
  )
}
