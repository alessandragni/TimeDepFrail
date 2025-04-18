% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check.formula_terms}
\alias{check.formula_terms}
\title{Check Correctness of Formula Terms}
\usage{
check.formula_terms(formula, data)
}
\arguments{
\item{formula}{Formula object specifying the relationship between the time-to-event, the covariates and the cluster variables.}

\item{data}{Dataset in which these variables can be found.}
}
\value{
An error if any condition is not satified.
}
\description{
The function controls that the terms composing the formula object provided in input to the main model are correct.
They must include:
\itemize{
\item response variable on the left hand side
\item covariates (numerical or categorical) on the right hand side
\item cluster variable (categorical) on the right hand side and specified by the function 'cluster()'
}

Moreover, it controls that the covariates are contained in the dataset provided.
}
\keyword{internal}
