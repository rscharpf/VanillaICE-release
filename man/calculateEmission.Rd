% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{calculateEmission}
\alias{calculateEmission}
\alias{calculateEmission,SummarizedExperiment-method}
\alias{calculateEmission,list-method}
\alias{calculateEmission,numeric-method}
\title{Calculate the emission probabilities for the 6-state HMM}
\usage{
calculateEmission(x, param = EmissionParam())
}
\arguments{
\item{x}{list of low-level data with two elements: a numeric vector of log
R ratios and a numeric vector of B allele frequencies}

\item{param}{parameters for the 6-state HMM}
}
\value{
A matrix of emission probabilities. Column correspond to
the HMM states and rows correspond to markers on the array (SNPs
and nonpolymorphic markers)
}
\description{
Given the data and an object containing parameters for the HMM,
this function computes emission probabilities.  This function is
not intended to be called by the user and is exported for internal
use by other BioC packages.
}
\seealso{
baumWelchUpdate
}

