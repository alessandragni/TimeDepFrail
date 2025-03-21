#' @title Extract Log-Likelihood for `AdPaik` Objects
#' @description Returns the log-likelihood of the fitted `AdPaik` model.
#' @param object An `AdPaik` model object.
#' @param ... Additional arguments (ignored).
#' @return A log-likelihood object with degrees of freedom (`df`).
#' @export
#' @importFrom stats logLik
logLik.AdPaik <- function(object, ...) {
  out <- object$Loglikelihood
  if (!is.null(object$NParameters)) attr(out, "df") <- object$NParameters
  else attr(out, "df") <- length(object$OptimalParameters)
  class(out) <- 'logLik'
  out
}

#-------------------------------------------------------------------------------

#' @title Extract Model Formula for `AdPaik`
#' @description Retrieves the formula used in the `AdPaik` model.
#' @param x An `AdPaik` model object.
#' @param ... Additional arguments (ignored).
#' @return The formula object.
#' @export
formula.AdPaik <- function(x, ...) {
  out <- x$formula
  out
}

#-------------------------------------------------------------------------------

#' @title Extract AIC for `AdPaik` Objects
#' @description Computes the AIC for an `AdPaik` model.
#' @param fit An `AdPaik` model object.
#' @param scale Not used (for compatibility).
#' @param k Penalty parameter (default is 2 for AIC).
#' @param ... Additional arguments (ignored).
#' @return A numeric vector with the number of parameters and AIC value.
#' @export
extractAIC.AdPaik <- function(fit, scale, k = 2, ...) {
  res <- logLik(fit)
  edf <- attr(res, "df")
  c(edf, -2 * res + k * edf)
}

#-------------------------------------------------------------------------------

#' @title Extract Number of Observations for `AdPaik`
#' @description Returns the number of observations used in the model.
#' @param object An `AdPaik` model object.
#' @param ... Additional arguments (ignored).
#' @return Integer: Number of observations.
#' @export
#' @importFrom stats nobs
nobs.AdPaik <- function(object, ...) {
  object$NObservations
}

#-------------------------------------------------------------------------------

#' @title 
#' Extracts the Coefficients for the 'Adapted Paik et Al.' Model
#'
#' @description
#' Extracts the optimal parameters  \eqn{\boldsymbol{\beta}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param object An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @details
#' The `coef.AdPaik` function extracts the coefficients from the 
#' `OptimalParameters` field in `object`.
#'
#' The function validates the structure of `object` and ensures compatibility 
#' with the expected model output. It throws an error if the object is malformed or 
#' inconsistent.
#'
#' @return A named list containing the coefficients.
#' 
#' @export
#'
#' @examples
#' # Example using the 'Academic Dropout' dataset
#' data(data_dropout)
#' 
#' # Define the formula and time axis for the model
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' eps <- 1e-10
#' categories_range_min <- c(-8, -2, eps, eps, eps)
#' categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)
#'
#'\donttest{
#' # Run the main model
#' result <- AdPaikModel(formula, data_dropout, time_axis,
#'                       categories_range_min, categories_range_max, TRUE)
#'
#' # Extract the coefficients
#' coef(result)
#' }

coef.AdPaik <- function (object, ...){
  
  # Check object structure
  check.result(object)
  
  # Extract information from input variables
  L <- n_intervals <- object$NIntervals
  R <- n_regressors <- object$NRegressors
  optimal_params <- object$OptimalParameters
  
  beta = optimal_params[(L+1):(L+R)]
  
  # Assign names to beta if regressors are provided
  if (!is.null(object$Regressors) && length(object$Regressors) == R) {
    names(beta) <- object$Regressors
  }
  
  class(beta) <- c("coef.AdPaik", class(beta))
  return (beta)
}


#-------------------------------------------------------------------------------
