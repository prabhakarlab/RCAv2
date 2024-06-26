% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateCellTypeFromProjection.R
\name{estimateCellTypeFromProjection}
\alias{estimateCellTypeFromProjection}
\title{Estimate the most likely cell type from the projection to the reference panel}
\usage{
estimateCellTypeFromProjection(
  rca.obj,
  confidence = NULL,
  ctRank = F,
  cSCompute = F
)
}
\arguments{
\item{rca.obj}{RCA object.}

\item{confidence}{a parameter indicating the difference between z-scores. If the difference is below this threshold, the cell type will be set to unknown (default NULL).}

\item{ctRank}{parameter indicating whether a relative rank coloring for each cell type shall be computed (default FALSE).}

\item{cSCompute}{parameter indicating wheter the confidence score should be computed for each cell (default FALSE).}
}
\value{
RCA object.
}
\description{
Estimate the most likely cell type from the projection to the reference panel
}
\examples{
RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
RCA.pbmcs <- dataClust(RCA.pbmcs)
RCA.pbmcs <- estimateCellTypeFromProjection(RCA.pbmcs)
}
