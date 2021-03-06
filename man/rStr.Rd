% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resHelpStrFns.R
\name{rStr}
\alias{rStr}
\title{Constructs an APA formatted string for a correlation, but from lm().}
\usage{
rStr(rVal, rdf, partial = F, ...)
}
\arguments{
\item{...}{Options.}

\item{param}{The parameter of interest. Can be "int" as a shortcut 
for "(Intercept)".}

\item{modObj}{Either an \code{lm} or \code{summary.lm} (faster) object.}
}
\description{
Constructs an APA formatted string for a correlation, but from lm().
}
