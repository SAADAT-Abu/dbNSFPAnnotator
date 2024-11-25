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
result <- annotate_variants(query_data, "path/to/dbNSFP.gz", chunk_size = 1000, workers = 6)
```
