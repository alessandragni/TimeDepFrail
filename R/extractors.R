

#-------------------------------------------------------------------------------
#' @title 
#' Extracts the Optimal Parameters of Each Category for the 'Adapted Paik et Al.' Model
#'
#' @description
#' Extracts the optimal parameters \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, 
#' \eqn{\nu}, \eqn{\boldsymbol{\gamma}} obtained with the 
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
#' @return A named list containing the categories of optimal parameters.
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
  
  phi = optimal_params[1:(L)]
  beta = optimal_params[(L+1):(L+R)]
  mu1 = optimal_params[(L+R+1)]
  nu = optimal_params[(L+R+2)]
  gamma = optimal_params[(L+R+3):(2*L+R+2)]
  
  
  # Assign names to beta if regressors are provided
  if (!is.null(object$Regressors) && length(object$Regressors) == R) {
    names(beta) <- object$Regressors
  }
  
  return_list = list('phi'=phi,
                     'beta'=beta,
                     'mu1'=mu1,
                     'nu'=nu,
                     'gamma'=gamma)
  class(return_list) <- c("coef.AdPaik", class(return_list))
  return (return_list)
}


#' @title 
#' Print Method for `coef.AdPaik` Objects
#'
#' @description
#' Displays the optimal parameter estimates extracted from an `AdPaik` model object in a structured format.
#'
#' @param x An object of class `coef.AdPaik`, typically returned by the `coef.AdPaik` function.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' This function prints the estimated optimal parameters for the 'Adapted Paik et al.' model.
#' It organizes the output into five categories: 
#' \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, \eqn{\nu}, and \eqn{\boldsymbol{\gamma}}.
#'
#' Each category is printed separately for clarity. If regressors are present, 
#' their corresponding \eqn{\beta} coefficients are displayed with names.
#'
#' @return Prints the coefficients to the console. No value is returned.
#' 
#' @export
print.coef.AdPaik <- function(x, ...){
  cat("Optimal Parameters:\n")
  cat("phi:\n")
  print(x$phi)
  cat("beta:\n")
  print(x$beta)
  cat("mu1:\n")
  print(x$mu1)
  cat("nu:\n")
  print(x$nu)
  cat("gamma:\n")
  print(x$gamma)
}


#-------------------------------------------------------------------------------
#' @title 
#' Extracts the Standard Errors Computed for Each Category for the 'Adapted Paik et Al.' Model
#'
#' @description
#' Extracts the standard errors for \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, 
#' \eqn{\nu}, \eqn{\boldsymbol{\gamma}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param object An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#'
#' @details
#' The `coefse` function extracts the standard errors for the estimated parameters from the 
#' `StandardErrorParameters` field in `object`. 
#'
#' The function validates the structure of `object` and ensures compatibility 
#' with the expected model output. It throws an error if the object is malformed or 
#' inconsistent.
#'
#' @return A named list containing the categories of the standard errors for the optimal parameters.
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
#' coefse(result)
#' }

coefse <- function(object){
  
  # Check object structure
  check.result(object)
  
  # Extract information from input variables
  L <- n_intervals <- object$NIntervals
  R <- n_regressors <- object$NRegressors
  optimal_params <- object$StandardErrorParameters
  
  phi = optimal_params[1:(L)]
  beta = optimal_params[(L+1):(L+R)]
  mu1 = optimal_params[(L+R+1)]
  nu = optimal_params[(L+R+2)]
  gamma = optimal_params[(L+R+3):(2*L+R+2)]
  
  
  # Assign names to beta if regressors are provided
  if (!is.null(object$Regressors) && length(object$Regressors) == R) {
    names(beta) <- object$Regressors
  }
  
  return_list = list('se.phi'=phi,
                     'se.beta'=beta,
                     'se.mu1'=mu1,
                     'se.nu'=nu,
                     'se.gamma'=gamma)
  class(return_list) <- c("coefse.AdPaik", class(return_list))
  return (return_list)
}


#-------------------------------------------------------------------------------
#' @title 
#' Print Method for `coefse.AdPaik` Objects
#'
#' @description
#' Displays the standard errors of the optimal parameter estimates extracted from an `AdPaik` model object.
#'
#' @param x An object of class `coefse.AdPaik`, typically returned by the `coefse` function.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' This function prints the standard errors for the estimated optimal parameters in the 'Adapted Paik et al.' model.
#' It organizes the output into five categories:
#' \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, \eqn{\nu}, and \eqn{\boldsymbol{\gamma}}.
#'
#' Each category is printed separately for clarity. If regressors are present, 
#' their corresponding \eqn{\beta} standard errors are displayed with names.
#'
#' @return Prints the standard errors to the console. No value is returned.
#' 
#' @export
print.coefse.AdPaik <- function(x, ...){
  cat("Optimal Parameters' standard errors:\n")
  cat("phi:\n")
  print(x$se.phi)
  cat("\n")
  cat("beta:\n")
  print(x$se.beta)
  cat("\n")
  cat("mu1:\n")
  print(x$se.mu1)
  cat("\n")
  cat("nu:\n")
  print(x$se.nu)
  cat("\n")
  cat("gamma:\n")
  print(x$se.gamma)
  cat("\n")
}




#-------------------------------------------------------------------------------
#' @title 
#' Extracts the Confidence Intervals Computed for Each Category for the 'Adapted Paik et Al.' Model
#'
#' @description
#' Extracts the confidence intervals for \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, 
#' \eqn{\nu}, \eqn{\boldsymbol{\gamma}} obtained with the 
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
#' The `confint.AdPaik` function extracts the standard errors for the estimated parameters from the 
#' `ParametersCI` field in `object`. 
#'
#' The function validates the structure of `object` and ensures compatibility 
#' with the expected model output. It throws an error if the object is malformed or 
#' inconsistent.
#'
#' @return A named list containing the categories of the standard errors for the optimal parameters.
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
#' confint(result)
#' }
confint.AdPaik <- function(object, parm = NULL, level = 0.95, ...){
  if(!is.null(parm))
      warning("Changing `parm` argument is not supported for this model. It will be ignored.")
  
  # Check object structure
  check.result(object)
  
  # Extract information from input variables
  L <- n_intervals <- object$NIntervals
  R <- n_regressors <- object$NRegressors
  
  optimal_params <- object$OptimalParameters
  se_optimal_params <- object$StandardErrorParameters
  confints = params_CI(optimal_params, se_optimal_params, level)
  
  optimal_params_left <- confints$ParamsCI_right
  
  phi_left = optimal_params_left[1:(L)]
  beta_left = optimal_params_left[(L+1):(L+R)]
  mu1_left = optimal_params_left[(L+R+1)]
  nu_left = optimal_params_left[(L+R+2)]
  gamma_left = optimal_params_left[(L+R+3):(2*L+R+2)]
  
  optimal_params_right <- confints$ParamsCI_right
  
  phi_right = optimal_params_right[1:(L)]
  beta_right = optimal_params_right[(L+1):(L+R)]
  mu1_right = optimal_params_right[(L+R+1)]
  nu_right = optimal_params_right[(L+R+2)]
  gamma_right = optimal_params_right[(L+R+3):(2*L+R+2)]
  
  
  # Assign names to beta if regressors are provided
  if (!is.null(object$Regressors) && length(object$Regressors) == R) {
    names(beta_left) <- object$Regressors
    names(beta_right) <- object$Regressors
  }
  
  phi = list("left" = phi_left, "right" = phi_right)
  beta = list("left" = beta_left, "right" = beta_right)
  mu1 = list("left" = mu1_left, "right" = mu1_right)
  nu = list("left" = nu_left, "right" = nu_right)
  gamma = list("left" = gamma_left, "right" = gamma_right)
  
  return_list = list('confint.phi'=phi,
                     'confint.beta'=beta,
                     'confint.mu1'=mu1,
                     'confint.nu'=nu,
                     'confint.gamma'=gamma)
  
  class(return_list) <- c("confint.AdPaik", class(return_list))
  return (return_list)
}
  
#-------------------------------------------------------------------------------
#' @title 
#' Print Method for `confint.AdPaik` Objects
#'
#' @description
#' Displays the confidence intervals for the optimal parameter estimates extracted from an `AdPaik` model object.
#'
#' @param x An object of class `confint.AdPaik`, typically returned by the `confint.AdPaik` function.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' This function prints the confidence intervals for the estimated optimal parameters in the 
#' 'Adapted Paik et al.' model. The output is organized into five categories:
#' \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, \eqn{\nu}, and \eqn{\boldsymbol{\gamma}}.
#'
#' Confidence intervals are displayed as two-column matrices (left and right bounds) where applicable.
#' Scalar values (\eqn{\mu_1} and \eqn{\nu}) are printed as individual numeric values.
#'
#' @return Prints the confidence intervals to the console. No value is returned.
#' 
#' @export
print.confint.AdPaik <- function(x, ...) {
  cat("Confidence Intervals for Optimal Parameters:\n")
  cat("phi:\n")
  print(do.call(cbind, x$confint.phi))
  cat("\n")
  cat("beta:\n")
  print(do.call(cbind, x$confint.beta))
  cat("\n")
  cat("mu1:\n")
  print(unlist(x$confint.mu1))
  cat("\n")
  cat("nu:\n")
  print(unlist(x$confint.nu))
  cat("\n")
  cat("gamma:\n")
  print(do.call(cbind, x$confint.gamma))
  cat("\n")
}






