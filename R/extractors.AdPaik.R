#' @title Extract Log-Likelihood for `AdPaik` Objects
#' @description Returns the log-likelihood of the fitted `AdPaik` model.
#' @param object An `AdPaik` model object.
#' @param ... Additional arguments (ignored).
#' @return A log-likelihood object with degrees of freedom (`df`).
#' @export
#' @method logLik AdPaik
logLik.AdPaik <- function(object, ...) {
  out <- object$Loglikelihood
  if (!is.null(result$NParameters)) attr(out, "df") <- result$NParameters
  else attr(out, "df") <- length(result$OptimalParameters)
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
#' @method formula AdPaik
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
#' @method extractAIC AdPaik
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
#' @method nobs AdPaik
nobs.AdPaik <- function(object, ...) {
  object$NObservations
}

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
