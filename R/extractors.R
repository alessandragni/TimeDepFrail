#-------------------------------------------------------------------------------
#' @title Extracts the optimal parameters of each cateogry for the 'Adapted Paik et al.' Model
#'
#' @description
#' Extracts the optimal parameters \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, 
#' \eqn{\nu}, \eqn{\boldsymbol{\gamma}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param result An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#'
#' @details
#' The `coef.AdPaik` function extracts the coefficients from the 
#' `OptimalParameters` field in the `result` object. 
#'
#' The function validates the structure of the `result` object and ensures compatibility 
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

coef.AdPaik <- function (result){
  
  # Check result structure
  check.result(result)
  
  # Extract information from input variables
  L <- n_intervals <- result$NIntervals
  R <- n_regressors <- result$NRegressors
  optimal_params <- result$OptimalParameters
  
  phi = optimal_params[1:(L)]
  beta = optimal_params[(L+1):(L+R)]
  mu1 = optimal_params[(L+R+1)]
  nu = optimal_params[(L+R+2)]
  gamma = optimal_params[(L+R+3):(2*L+R)]
  
  
  # Assign names to beta if regressors are provided
  if (!is.null(result$Regressors) && length(result$Regressors) == R) {
    names(beta) <- result$Regressors
  }
  
  return_list = list('phi'=phi,
                     'beta'=beta,
                     'mu1'=mu1,
                     'nu'=nu,
                     'gamma'=gamma)
  return (return_list)
}


#-------------------------------------------------------------------------------
#' @title Extracts the standard errors computed for each cateogry for the 'Adapted Paik et al.' Model
#'
#' @description
#' Extracts the standard errors for \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, 
#' \eqn{\nu}, \eqn{\boldsymbol{\gamma}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param result An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#'
#' @details
#' The `coefse.AdPaik` function extracts the standard errors for the estimated parameters from the 
#' `StandardErrorParameters` field in the `result` object. 
#'
#' The function validates the structure of the `result` object and ensures compatibility 
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

coefse <- function(result){
  
  # Check result structure
  check.result(result)
  
  # Extract information from input variables
  L <- n_intervals <- result$NIntervals
  R <- n_regressors <- result$NRegressors
  optimal_params <- result$StandardErrorParameters
  
  phi = optimal_params[1:(L)]
  beta = optimal_params[(L+1):(L+R)]
  mu1 = optimal_params[(L+R+1)]
  nu = optimal_params[(L+R+2)]
  gamma = optimal_params[(L+R+3):(2*L+R)]
  
  
  # Assign names to beta if regressors are provided
  if (!is.null(result$Regressors) && length(result$Regressors) == R) {
    names(beta) <- result$Regressors
  }
  
  return_list = list('se.phi'=phi,
                     'se.beta'=beta,
                     'se.mu1'=mu1,
                     'se.nu'=nu,
                     'se.gamma'=gamma)
  return (return_list)
}



#-------------------------------------------------------------------------------
#' @title Extracts the confidence intervals computed for each cateogry for the 'Adapted Paik et al.' Model
#'
#' @description
#' Extracts the confidence intervals for \eqn{\boldsymbol{\phi}}, \eqn{\boldsymbol{\beta}}, \eqn{\mu_1}, 
#' \eqn{\nu}, \eqn{\boldsymbol{\gamma}} obtained with the 
#' time-dependent frailty model proposed in the 'Adapted Paik et al.' framework.
#'
#' @param result An S3 object of class `AdPaik`, returned by the main model function 
#' (`AdPaikModel`). This object contains all the optimal parameter estimates.
#'
#' @details
#' The `confint.AdPaik` function extracts the standard errors for the estimated parameters from the 
#' `ParametersCI` field in the `result` object. 
#'
#' The function validates the structure of the `result` object and ensures compatibility 
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

confint.AdPaik <- function(result){
  
  # Check result structure
  check.result(result)
  
  # Extract information from input variables
  L <- n_intervals <- result$NIntervals
  R <- n_regressors <- result$NRegressors
  optimal_params_left <- result$ParametersCI$ParamsCI_left
  
  phi_left = optimal_params_left[1:(L)]
  beta_left = optimal_params_left[(L+1):(L+R)]
  mu1_left = optimal_params_left[(L+R+1)]
  nu_left = optimal_params_left[(L+R+2)]
  gamma_left = optimal_params_left[(L+R+3):(2*L+R)]
  
  optimal_params_right <- result$ParametersCI$ParamsCI_right
  
  phi_right = optimal_params_right[1:(L)]
  beta_right = optimal_params_right[(L+1):(L+R)]
  mu1_right = optimal_params_right[(L+R+1)]
  nu_right = optimal_params_right[(L+R+2)]
  gamma_right = optimal_params_right[(L+R+3):(2*L+R)]
  
  
  # Assign names to beta if regressors are provided
  if (!is.null(result$Regressors) && length(result$Regressors) == R) {
    names(beta_left) <- result$Regressors
    names(beta_right) <- result$Regressors
  }
  
  phi = list("left" = phi_left, "right" = phi_right)
  beta = list("left" = beta_left, "right" = beta_right)
  mu1 = list("left" = mu1_left, "right" = mu1_right)
  nu = list("left" = nu_left, "right" = nu_right)
  gamma = list("left" = gamma_left, "right" = gamma_right)
  
  class(phi) <- "phi"
  class(beta) <- "beta"
  class(mu1) <- "mu1"
  class(nu) <- "nu"
  class(gamma) <- "gamma"
  
  return_list = list('confint.phi'=phi,
                     'confint.beta'=beta,
                     'confint.mu1'=mu1,
                     'confint.nu'=nu,
                     'confint.gamma'=gamma)
  return (return_list)
}




