% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isos_err_distrib.R
\name{isos_err_distrib}
\alias{isos_err_distrib}
\title{Sample in a distribution of values for which an error value in associated.}
\usage{
isos_err_distrib(vc_values, vc_errors, sample_size)
}
\arguments{
\item{vc_values}{List of values to be sampled.}

\item{vc_errors}{List of error values to be sampled.}

\item{sample_size}{size of the samples to be drawn.}
}
\value{
A vector of size \code{sample_size}.
}
\description{
Sample in a distribution of values for which an error value in associated.
}
\examples{
isos_err_distrib(1:10, rep(2,10), 100)
}
