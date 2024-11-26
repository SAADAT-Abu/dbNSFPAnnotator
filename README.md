# dbNSFPAnnotator

## Overview
`dbNSFPAnnotator` is an R package for annotating variants using the dbNSFP database with parallel processing capabilities.

## Installation
To install this package from GitHub:
```R
# Install devtools if not already installed
install.packages("devtools")

# Install dbNSFPAnnotator from GitHub
devtools::install_github("SAADAT-Abu/dbNSFPAnnotator")
```

## Usage

```R
library(dbNSFPAnnotator)

# Example query data
query_data <- data.frame(
  chr = c("1", "2"),
  start = c(12345, 67890),
  end = c(12345, 67890),
  ref = c("A", "G"),
  alt = c("T", "C")
)

# Annotate variants using a bgzipped dbNSFP file

annotate_variants(
  query_data,          # Data frame with variant information
  dbnsfp_file,         # Path to the bgzipped dbNSFP file
  columns = NULL,  # Vector of dbNSFP columns to include (default: all)
  chunk_size = 1000,   # Number of variants to process per chunk
  workers = 6          # Number of cores for parallel processing
)

```

## Example

```R

dbnsfp_file <- "path/to/dbNSFP.gz"

# List all available columns
columns <- dbNSFP_columns(dbnsfp_file)
print(columns)

# Prepare query data
query_data <- data.frame(
  chr = c("1", "2"),
  start = c(12345, 67890),
  end = c(12345, 67890),
  ref = c("A", "G"),
  alt = c("T", "C")
)

# Annotate using specific columns
selected_columns <- c("chr", "pos(1-based)", "SIFT_score", "Polyphen2_HDIV_score")
result <- annotate_variants(query_data, dbnsfp_file, selected_columns = selected_columns)

print(result)

# Annotate using all columns
result_all <- annotate_variants(query_data, dbnsfp_file)
print(result_all)

```


