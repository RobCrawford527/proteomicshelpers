% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_reps.R
\name{combine_reps}
\alias{combine_reps}
\title{Calculate Mean Across Replicates}
\usage{
combine_reps(data, samples = NULL, sam_col, rep_col, val_col, min_reps = 1)
}
\arguments{
\item{data}{A data frame containing multiple replicates}

\item{samples}{List of samples to calculate means for (defaults to all)}

\item{sam_col}{Column indicating the samples}

\item{rep_col}{Column indicating the replicate}

\item{val_col}{Column containing the values}

\item{min_reps}{Minimum number of replicates a protein must be present
in to have mean calculated (default 1)}
}
\value{
A data frame in the same format as the input, containing the
mean values
}
\description{
Calculate Mean Across Replicates
}
