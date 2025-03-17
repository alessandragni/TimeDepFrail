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
#' Extract Processed Dataset from an AdPaik Object
#' @description 
#' Retrieves the dataset used in the 'Adapted Paik et al.' model, applying any 
#' missing value handling specified in the model object.
#' @param object An `AdPaik` model object.
#' @return The processed dataset used in the model, or `NULL` if no dataset is found.
#' @export
#' @importFrom nlme getData
getData.AdPaik <- function(object) {
  data <- object$dataset
  if (is.null(data)) return(NULL)
  
  # Handle missing values if `na.action` is present
  naAct <- object[["na.action"]]
  if (!is.null(naAct)) {
    if (inherits(naAct, "omit")) 
      data <- data[-naAct, ]
    else if (inherits(naAct, "exclude")) 
      data <- data
    else 
      data <- eval(naAct)(data)
  }
  
  return(data)
}

#-------------------------------------------------------------------------------
