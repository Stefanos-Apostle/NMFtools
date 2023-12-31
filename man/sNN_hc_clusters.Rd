% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF_utils.R
\name{sNN_hc_clusters}
\alias{sNN_hc_clusters}
\title{sNN HClust Clusters}
\usage{
sNN_hc_clusters(
  sNN_clusters,
  trait_dists,
  clusters_to_calc = "All",
  method = "euclidean"
)
}
\value{
Hclust object ordering the sNN clusters
}
\description{
Function to calc hclust of sNN clusters, as a way to compare shared nearest neighbors cluters.
}
\examples{
sNN_hc_clusters(...)
}
