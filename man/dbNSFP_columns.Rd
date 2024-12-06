% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbNSFP_columns.R
\name{dbNSFP_columns}
\alias{dbNSFP_columns}
\title{List Available Columns in dbNSFP}
\usage{
dbNSFP_columns(dbnsfp_file)
}
\arguments{
\item{dbnsfp_file}{Path to the bgzipped dbNSFP file.}
}
\value{
A character vector of column names if the file is valid and has a proper header.
}
\description{
This function lists all available column names in a dbNSFP file.
}
