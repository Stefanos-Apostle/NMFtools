% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF_utils.R
\name{goi_manhattan}
\alias{goi_manhattan}
\title{Gene of interest manhattan plot}
\usage{
goi_manhattan(top_bmi_genes, NMF, factor)
}
\arguments{
\item{top_bmi_genes}{Genes to plot (Need to change this object name)}

\item{NMF}{NMF object}

\item{factor}{NMF Factor to plot}
}
\value{
GGplot of contribution manhattans by genes for a factor
}
\description{
Function to plot manhattans of genes of interest for a specific factor, faceted by nearest gene
}
\examples{
goi_manhattan(top_bmi_genes = c("ADIPOQ", "FTO", "LEP"), NMF = NMF_object)
}
