% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF_utils.R
\name{fix_pval_float}
\alias{fix_pval_float}
\title{Fix Pval Float}
\usage{
fix_pval_float(enrichment_df)
}
\arguments{
\item{enrichment_df}{Output of enrichment_calc}
}
\value{
Updated enrichment_df
}
\description{
Function to replace pval = 0 with smallest r value
}
\examples{
enrichment_df <- fix_pval_float(enirchment_df)
}
