% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_faust_pop.R
\name{.get_faust_pop}
\alias{.get_faust_pop}
\title{Filter an expression matrix to return only cells matching one or more FAUST annotation combinations}
\usage{
.get_faust_pop(ex, sample, project_path, pop)
}
\arguments{
\item{ex}{\code{matrix}. Matrix containing marker expression values for \code{sample}.}

\item{sample}{\code{character}. Name of sample, as found in directory
<project_path>/faustData/sampleData/.}

\item{project_path}{character. Path to directory containing FAUST output.}

\item{pop}{\code{list} or \code{named character vector}. If a \code{character vector},
then all cells matching the set of marker levels are returned. If a \code{list}, then each
element must be a \code{character vector}, and then
all cells matching either of these specified sets of marker levels are returned.}
}
\value{
Numeric matrix.
}
\description{
Filter an expression matrix to return only cells matching one or more FAUST annotation combinations
}