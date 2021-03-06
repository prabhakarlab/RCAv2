% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createRCAObjectFrom10X.R
\name{createRCAObjectFrom10X}
\alias{createRCAObjectFrom10X}
\title{Create RCA Object from 10X data}
\usage{
createRCAObjectFrom10X(dataDir, min.barcode.umi = 100)
}
\arguments{
\item{dataDir}{Directory containing 10X data}

\item{min.barcode.umi}{Minimum UMIs needed for a barcode to be considered (default 100).}
}
\value{
RCA object.
}
\description{
Create RCA Object from 10X data
}
\examples{
\dontrun{
RCA.pbmcs <- createRCAObjectFrom10X("directory_containing_10X_files",min.barcode.umi=50)
print(RCA.pbmcs)
}
}
