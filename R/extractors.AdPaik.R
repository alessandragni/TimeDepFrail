
#-------------------------------------------------------------------------------

#' @title Extract Number of Observations for `AdPaik`
#' @description Returns the number of observations used in the model.
#' @param object An `AdPaik` model object.
#' @param ... Additional arguments (ignored).
#' @return Integer: Number of observations.
#' @importFrom stats nobs
#' @export
nobs.AdPaik <- function(object, ...) {
  object$NObservations
}

#-------------------------------------------------------------------------------

#' @title Extract Log-Likelihood for `AdPaik` Objects
#' @description Returns the log-likelihood of the fitted `AdPaik` model.
#' @param object An `AdPaik` model object.
#' @param ... Additional arguments (ignored).
#' @return A log-likelihood object with degrees of freedom (`df`).
#' @importFrom stats logLik
#' @export
logLik.AdPaik <- function(object, ...) {
  out <- object$Loglikelihood
  if (!is.null(object$NParameters)) attr(out, "df") <- object$NParameters
  else attr(out, "df") <- length(object$OptimalParameters)
  class(out) <- 'logLik'
  out
}


#-------------------------------------------------------------------------------

#' @title Extract AIC for `AdPaik` Objects
#' @description Computes the AIC for an `AdPaik` model.
#' @param fit An `AdPaik` model object.
#' @param scale Changing it is not supported for this model. It will be ignored.
#' @param k Penalty parameter (default is 2 for AIC).
#' @param ... Additional arguments (ignored).
#' @return A numeric vector with the number of parameters and AIC value.
#' @importFrom stats extractAIC
#' @export
extractAIC.AdPaik <- function(fit, scale = NULL, k = 2, ...) {
  res <- logLik(fit)
  edf <- attr(res, "df")
  c(edf, -2 * as.numeric(res) + k * edf)
}

#-------------------------------------------------------------------------------

#' @title 
#' Extract the Coefficients for the 'Adapted Paik et Al.' Model
#'
#' @description
#' Extracts the optimal \eqn{\boldsymbol{\beta}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param object An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#' @param ... Additional arguments (ignored).
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
  
  return (beta)
}




#-------------------------------------------------------------------------------
#' @title 
#' Extracts the Confidence Intervals for the Coefficients for the 'Adapted Paik et Al.' Model
#'
#' @description
#' Extracts the confidence intervals for \eqn{\boldsymbol{\beta}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param object An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#' @param parm A specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. 
#' Defaults to NULL, and all parameters are considered. Changing it is not supported for this model. It will be ignored.
#' @param level The confidence level required. Defaults to 0.95.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @details
#' The `confint.AdPaik` function extracts the standard errors for the beta coefficients from the 
#' `ParametersCI` field in `object`. 
#'
#' The function validates the structure of `object` and ensures compatibility 
#' with the expected model output. It throws an error if the object is malformed or 
#' inconsistent.
#'
#' @return A named list containing the categories of the standard errors for the optimal parameters.
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
#' confint(result)
#' }
confint.AdPaik <- function(object, parm = NULL, level = 0.95, ...) {
  if (!is.null(parm)) 
    warning("Changing `parm` argument is not supported for this model. It will be ignored.")
  
  # Validate object structure
  check.result(object)
  
  # Extract information
  L <- object$NIntervals
  R <- object$NRegressors
  optimal_params <- object$OptimalParameters
  se_optimal_params <- object$StandardErrorParameters
  
  # Compute confidence intervals
  confints <- params_CI(optimal_params, se_optimal_params, level)
  beta_left <- confints$ParamsCI_left[(L + 1):(L + R)]
  beta_right <- confints$ParamsCI_right[(L + 1):(L + R)]
  
  # Assign names if regressors are available
  if (!is.null(object$Regressors) && length(object$Regressors) == R) {
    names(beta_left) <- object$Regressors
    names(beta_right) <- object$Regressors
  }
  
  # Define confidence level percent labels
  a <- (1 - level) / 2
  pct <- format_perc(c(a, 1 - a), 3)
  
  # Create output matrix
  ci <- matrix(NA_real_, nrow = R, ncol = 2, 
               dimnames = list(names(beta_left), pct))
  
  ci[, 1] <- beta_left
  ci[, 2] <- beta_right
  
  return(ci)
}




#-------------------------------------------------------------------------------

#' Plots Related to the the 'Adapted Paik et Al.' Model
#'
#' @param x An object of class 'AdPaik'.
#' @param which A numeric vector indicating which plots to display. 
#'             Choices: 1 = Baseline Hazard, 
#'                      2 = Posterior Frailty Estimate.
#' @param captions A character vector with captions for each plot.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @return No return value. This function generates plots.
#' 
#' @examples
#' # Import data
#' data(data_dropout)
#' 
#' # Define the variables needed for the model execution
#' eps_paik <- 1e-10
#' categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
#' categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#'
#' # Call the main model function
#' \donttest{
#' result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#'
#' plot(result)
#' }  
#' 
#' @export

plot.AdPaik <- function(x, which = c(1, 2), 
                        captions = c("Plot 1: Baseline Hazard", 
                                     "Plot 2: Posterior Frailty Estimate"), ...) {
  
  if (!inherits(x, "AdPaik")) 
    stop("This function is only applicable to 'AdPaik' objects")
  
  if (!is.numeric(which) || any(which < 1) || any(which > 2)) 
    stop("'which' must be in 1:2")
  
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  
  # Plot 1:
  if (show[1]) {
    plot_bas_hazard(x)
    grDevices::dev.flush()
  }
  
  # Plot 2: 
  if (show[2]) {
    plot_post_frailty_est(x)
    grDevices::dev.flush()
  }
  
  
  invisible()
}



