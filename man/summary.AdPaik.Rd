% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.AdPaik}
\alias{summary.AdPaik}
\title{Summary Method for AdPaik Objects}
\usage{
\method{summary}{AdPaik}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"AdPaik"}, returned by the main model function.}

\item{...}{Additional arguments (currently unused).}
}
\value{
\code{summary.AdPaik()} returns an object of class \code{"summary.AdPaik"}, containing structured model summary information.
}
\description{
\code{summary} method for objects of class \code{"AdPaik"}. Provides a structured summary of the results from the Adapted Paik et al.'s Time-Dependent Shared Frailty Model.
}
\details{
This function extracts relevant model outputs, such as estimated parameters, standard errors,
and the convergence status of the algorithm. The summary is returned as an object of class \code{"summary.AdPaik"},
which can be printed using \code{print()}.
}
\examples{
\dontrun{
data(data_dropout)
eps_paik <- 1e-10
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
formula <- time_to_event ~ Gender + CFUP + cluster(group)

# Fit the model
result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)

# Generate and print summary
summary_result <- summary(result)
print(summary_result)  # or simply `summary(result)`
}
}
\seealso{
\code{\link{summary}}, \code{\link{print}}
}
