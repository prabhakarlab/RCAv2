% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataClust.R
\name{dataClust}
\alias{dataClust}
\title{Generate cell clusters using hierarchical clustering and dynamic tree cutting.}
\usage{
dataClust(rca.obj, deepSplitValues = 1, minClustSize = 5, corMeth = "pearson")
}
\arguments{
\item{rca.obj}{RCA object.}

\item{deepSplitValues}{integer vector indicating how deep the dendrogram should be cut. Values can range from 0 to 4 (default 1).}

\item{minClustSize}{integer value indicating the minimum size of the resulting clusters (default 5).}

\item{corMeth}{Correlation method used to compute the distance matrix of the projection (pearson (default), spearman, kendal).}
}
\value{
RCA object.
}
\description{
Generate cell clusters using hierarchical clustering and dynamic tree cutting.
}
\examples{
RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
RCA.pbmcs <- dataClust(RCA.pbmcs)

}
