% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/params_se.R
\name{params_CI}
\alias{params_CI}
\title{Confidence interval for the optimal estimated parameters}
\usage{
params_CI(optimal_params, se_params)
}
\arguments{
\item{optimal_params}{Numerical vector of optimal estimated parameters. Its length is equal to the number of model parameters.}

\item{se_params}{Numerical vector containing the standard error associated to each estimated parameter.}
}
\value{
A S3  object of class 'ParametersCI', composed of two numerical vector of length equal to the number of model parameters:
\itemize{
\item ParamsCI_left: left confidence interval for each parameter
\item ParamsCI_right: right confidence interval for each parameter
}
}
\description{
The function provides the confidence interval for each estimated parameter, using the standard error
computed through another method and provided as second argument to the current function.
}
