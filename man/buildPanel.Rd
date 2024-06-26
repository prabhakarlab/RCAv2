% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildReferencePanel.R
\name{buildPanel}
\alias{buildPanel}
\title{Internal function to build a reference panel from your bulk RNA sequencing data and a list of marker genes
Bulk RNA sequencing data must contain gene symbols - we recommend conversion prior to use this function, e.g. using Biomart.}
\usage{
buildPanel(bulk.rna.data, de.genes, species = "HUMAN")
}
\arguments{
\item{bulk.rna.data}{Bulk RNA sequencing data - ideally TPM normalized}

\item{de.genes}{list of marker genes}

\item{species}{"HUMAN" for human data, "MOUSE" for mouse data. Default is "HUMAN".}
}
\value{
reference panel boject
}
\description{
Internal function to build a reference panel from your bulk RNA sequencing data and a list of marker genes
Bulk RNA sequencing data must contain gene symbols - we recommend conversion prior to use this function, e.g. using Biomart.
}
