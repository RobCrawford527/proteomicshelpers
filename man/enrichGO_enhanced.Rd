% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichGO_enhanced.R
\name{enrichGO_enhanced}
\alias{enrichGO_enhanced}
\title{Perform Gene Ontology Analysis}
\usage{
enrichGO_enhanced(
  genes,
  OrgDb,
  keyType = "ORF",
  ont = "all",
  reference,
  simplify = TRUE,
  evaluate = TRUE,
  convert = TRUE,
  toType = "GENENAME"
)
}
\arguments{
\item{genes}{The vector of gene names to be analysed}

\item{OrgDb}{Organism database to use}

\item{keyType}{Starting format of the names (default UNIPROT), as
defined by clusterProfiler}

\item{ont}{The ontologies to analyse (BP, CC and/or MF; defaults to "all")}

\item{reference}{Data frame containing the background list in multiple
formats}

\item{simplify}{Whether to simplify the GO result}

\item{evaluate}{Whether to evaluate the GO result}

\item{convert}{Whether to convert the names in the GO result}

\item{toType}{Format to convert the names to (if convert = TRUE)}
}
\value{
A data frame containing the enriched GO terms
}
\description{
Perform Gene Ontology Analysis
}
