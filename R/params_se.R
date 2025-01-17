#' @title
#' Standard error of the parameters
#'
#' @description
#' Function for computing the standard error of each optimal parameter, estimated through the
#' constraint multi-dimensional optimization.
#' The procedure for the computation is based on the numerical approximation of the second derivative
#' of the log-likelihood function, by the 'centered finite difference scheme' with an accuracy of the second order.
#'
#' @details
#' The standrd error of each parameter is computed as the inverse of the square root of the 'Information matrix', that in turn
#' is computed as the opposite of the 'Hessian matrix'. Only its diagonal is built and its elements are separatey
#' evaluated through a numerical approximation of the second derivative of the log-likelihood function.
#'
#' The function requires the optimal parameter vector and other parameters-related variables, to check:
#' - the right numerosity of the parameter vector
#' - the correct range existence of each parameter (i.e. each parameter lies in its range).
#'
#' @param optimal_params Numerical vector of optimal parameters. Its length (i.e. number of parameters) is equal to \eqn{n_p}.
#' @param params_range_min Numerical vector of length equal to \eqn{n_p}, containing the minimum range of each parameter.
#' @param params_range_max Numerical vector of length equal to \eqn{n_p}, containing the maximum range of each parameter.
#' @param dataset Dataset containing the value of the regressors for all individuals in the study.
#' @param centre vector containing the group membership of each individual and that induces the clustering subdivision.
#' @param time_axis Temporal domain.
#' Its number of intervals corresponds to the length of the time-domain minus 1
#' @param dropout_matrix Binary matrix of dimension (n_individuals, n_intervals).
#' The sum of the elements of each row must be (1), if the associated individual failed in a precise interval, and (0) if the individual
#' did not fail in the @time-axis.
#' Therefore, if an individual failed in the time-domain, the interval in which he failed will have value (1) and the others (0).
#' @param e_matrix Matrix of dimension (n_individuals, n_intervals) where each element contains the resolution of the temporal
#' integral for that individual in that interval, thorugh the 'e_time_fun' function.
#' @param h_dd Discretization step for the numerical approximation of the second derivative fo the loglikelihood function.
#'
#' @return Vector of parameter standard error, of length equal to the number of model parameters.

params_se.AdPaik <- function(optimal_params, params_range_min, params_range_max,
                     dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd){
  
  # Extract information from input variables
  n_params <- length(optimal_params)
  n_intervals <- length(time_axis) - 1
  n_regressors <- dim(dataset)[2]
  
  # Initialize parameters standard error vector
  se <- rep(0, n_params)
  
  # Initialize information and hessian element
  information_element <- hessian_element <- 0
  
  # Compute parameters se
  for(p in 1:n_params){
    # Store current parameter value and saved its updated value
    value <- optimal_params[p]
    value_plus_h <- value + h_dd
    value_minus_h <- value - h_dd
    if(value_plus_h > params_range_max[p])
      values_plus_h <- params_range_max[p]
    else if(value_minus_h < params_range_min[p])
      value_minus_h <- params_range_min[p]
    
    # Store original optimal parameters
    params_plus  <- optimal_params
    params_minus <- optimal_params
    
    # Update current parameters value
    params_plus[p] <- value_plus_h
    params_minus[p] <- value_minus_h
    
    # Compute log-likelihood function in new values and current value
    ll_eval <- ll_AdPaik_eval(optimal_params, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_eval_plus <- ll_AdPaik_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_eval_minus <- ll_AdPaik_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    
    # Approximate the second derivative of the log-likelihood function
    hessian_element <- (ll_eval_plus + ll_eval_minus - 2*ll_eval)/(h_dd * h_dd)
    
    if((hessian_element == Inf) || (hessian_element == -Inf)){
      #information_element <- hessian_element
      se[p] <- 1e-4
    }
    else{
      # Compute the information element from the hessian
      information_element <- - hessian_element
      
      # Compute standard error of the parameter
      se[p] <- 1/sqrt(information_element)
    }
  }
  
  # Return the entire vector of parameter standard error
  return (se)
}

#-------------------------------------------------------------------------------
#' @title
#' Confidence interval for the optimal estimated parameters
#'
#' @description
#' The function provides the confidence interval for each estimated parameter, using the standard error
#' computed through another method and provided as second argument to the current function.
#'
#'
#' @param optimal_params Numerical vector of optimal estimated parameters. Its length is equal to the number of model parameters.
#' @param se_params Numerical vector containing the standard error associated to each estimated parameter.
#' @param level A numeric value representing the confidence level.
#'
#' @return A S3  object of class 'ParametersCI', composed of two numerical vector of length equal to the number of model parameters:
#' - ParamsCI_left: left confidence interval for each parameter
#' - ParamsCI_right: right confidence interval for each parameter

params_CI <- function(optimal_params, se_params, level){
  
  # Check both input variables have the same dimension
  if(length(optimal_params) != length(se_params))
    stop("'optimal_params' and 'se_params' have different length.")
  
  # Check level
  if(level>1 || level<0)
    stop("'level' should be between 0 and 1")
  
  # Extract information from input variables
  n_params <- length(optimal_params)
  
  # Compute critical z-score for the given confidence level
  alpha <- 1 - level
  z_critical <- stats::qnorm(1 - alpha / 2)  # Two-tailed
  
  # Confidence interval for the parameters
  params_CI_left <- params_CI_right <- rep(0, n_params)
  for(p in 1:n_params){
    params_CI_left[p] <- optimal_params[p] - z_critical * se_params[p]
    params_CI_right[p] <- optimal_params[p] + z_critical * se_params[p]
  }
  
  params_CI <- list("ParamsCI_left" = params_CI_left,
                    "ParamsCI_right" = params_CI_right)
  class(params_CI) <- "ParametersCI"
  return (params_CI)
}
#-------------------------------------------------------------------------------
