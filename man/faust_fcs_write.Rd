% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_write.R
\name{faust_fcs_write}
\alias{faust_fcs_write}
\title{Save FAUST subset as an FCS file}
\usage{
faust_fcs_write(
  project_path,
  dir_save = NULL,
  pop,
  fr_source = NULL,
  sample = NULL,
  trans_fn = NULL,
  trans_chnl = NULL
)
}
\arguments{
\item{project_path}{character. FAUST project directory.}

\item{dir_save}{character.
Directory to save to.
If \code{NULL} (the default), then FCS files
are saved to
\code{file.path(project_path, "faustData", "fcsData", <pop_defn>)},
where \code{<pop_defn>} is a concatenation of the population
definition, e.g list(c("CD3" = 1, "CD4" = 2)) becomes "CD3~1~CD4~2~".
If \code{character}, then the FCS files are saved directly to this directory.
As stated before, default is \code{NULL}.}

\item{pop}{\code{list}, \code{named character vector}
or unnamed character vector of length one.
If an unnamed character vector,
then all cells matching that full FAUST annotation will be returned
(e.g. "CD3+CD8+CD45RA+IFNg+IL2-TNF+" or, equivalently (depending
on the number of sub-populations FAUST detects)
"CD3~2~2~CD8~2~2~CD45RA~2~2~IFNg~2~2~IL2~1~2~TNF~2~2~")).
If a named \code{character vector},
then all cells matching the set of marker levels are returned.
(names are markers, and elements are levels, e.g. c("CD3" = 2, "CD8" = 2)).
If a \code{list}, then each
element must be a named \code{character vector}, and then
all cells matching either of
these specified sets of marker levels are returned.}

\item{fr_source}{GatingSet or character vector. If a \code{GatingSet}, then
the flowFrames within the GatingSet are used to create the output fcs files.
If a character vector,
then must specify a directory containing FCS files. These FCS
files are used to create the output fcs files.}

\item{sample}{integer vector or character vector.
If an integer vector, then it
specifies the indices of the samples for which to save output. If character,
then it specifies the names of the sames to save output for. If \code{NULL},
then the output for every sample is saved. Default is \code{NULL}.}

\item{trans_fn}{function. If supplied, this function
is applied to the expression data.
Useful for back-transformation.
If \code{NULL}, then no transformation is applied. Default is \code{NULL}.}

\item{trans_chnl}{character vector.
If specified, \code{trans_fn} is applied to only these channels.
If \code{NULL} and if \code{trans_fn} is not \code{NULL},
then \code{trans_fn} is applied to entire
expression matrix. Default is \code{NULL}.}
}
\value{
\code{invisible(TRUE)}. Side effect is the saved FCS file.
}
\description{
Save a specified FAUST-identified population to
\code{project_path/faustData/fcsData}
as an FCS file for all samples gated or just a specified subset.
}
\examples{
faust_fcs_write(
  project_path = "", pop = list("CD3" = 2),
  gs = gs, sample = 1
)

}
