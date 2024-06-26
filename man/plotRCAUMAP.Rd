% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotRCAUMAP.R
\name{plotRCAUMAP}
\alias{plotRCAUMAP}
\title{Plot umap of projection to the RCA panel}
\usage{
plotRCAUMAP(
  rca.obj,
  cellPropertyList = NULL,
  filename = "RCA_UMAP.pdf",
  fontsize = 10
)
}
\arguments{
\item{rca.obj}{RCA object}

\item{cellPropertyList}{list of cell properties to plot}

\item{filename}{file name of saved umap}

\item{fontsize}{Size of the font used for plotting}
}
\value{
list of UMAP plots
}
\description{
The presence of cell type estimates, relative ranks and confindence scores are detected automatically and are plotted accordingly.
}
\examples{
\dontrun{
RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
RCA.pbmcs <- computeUMAP(RCA.pbmcs)
plotRCAUMAP(RCA.pbmcs)
}

}
