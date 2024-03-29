% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_init_exp.R
\name{create_init_exp}
\alias{create_init_exp}
\title{Single-cell SummarizedExperiment object}
\usage{
create_init_exp(counts, coldata)
}
\arguments{
\item{counts}{Named numeric matrix with features in rows and samples in columns.}

\item{coldata}{A data frame containing the annotations of the samples.}
}
\value{
SummarizedExperiment with provided data
}
\description{
Creates a SummarizedExperiment object necessary to build a multi-view representation.
}
\details{
This function is the first step for a multicellular factor analysis.
It collects in a single object the pseudobulk counts of a single cell experiment
and its annotations.
}
\examples{
inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
load(file.path(inputs_dir, "testpbcounts.rda"))
load(file.path(inputs_dir, "testcoldata.rda"))
pb_obj <- create_init_exp(counts = testpbcounts,  coldata = testcoldata)
}
