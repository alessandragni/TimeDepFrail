#' @title
#' Baseline Hazard Step-Function
#'
#' @description
#' The method computes the baseline hazard step-function in each interval of the time-domain, using the estimated parameters
#' \eqn{\phi_k, \forall k}
#'
#' @param object S3 object of class 'AdPaik' returned by the main model output, that contains all the information for the computation
#' of the frailty standard deviation.
#'
#' @return Numerical vector of length equal to the number of intervals of the time-domain, with the value of the baseline hazard step-function.
#' 
#' @export 
#' 
#' @examples
#' # Consider the 'Academic Dropout dataset'
#' data(data_dropout)
#'
#' # Define the variables needed for the model execution
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' eps <- 1e-10
#' categories_range_min <- c(-8, -2, eps, eps, eps)
#' categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)
#'
#' \donttest{
#' # Call the main model
#' result <- AdPaikModel(formula, data_dropout, time_axis,
#'                       categories_range_min, categories_range_max)
#'
#' bas_hazard(result)
#' }
bas_hazard <- function(object){
  
  # Return the baseline hazard step-function
  return(object$BaselineHazard)
}

#_______________________________________________________________________________________________________________________

#' @title
#' Internal Function for the Baseline Hazard Step-Function
#'
#' @description
#' The method computes the baseline hazard step-function in each interval of the time-domain, using the estimated parameters
#' \eqn{\phi_k, \forall k}
#'
#'
#' @param optimal_params Numerical vector of length equal to the number of model parameters, containing the optimal estimated parameters.
#' @param time_axis Numerical vector of temporal domain.
#'
#' @return Numerical vector of length equal to the number of intervals of the time-domain, with the value of the baseline hazard step-function.
#'
#' @keywords internal
#' 
bas_hazard_internal <- function(optimal_params, time_axis){
  # Extract information from input variables
  L <- n_intervals <- length(time_axis) - 1
  eps <- 1e-2

  # Extract baseline hazard parameters from input parameters
  phi <- optimal_params[1:L]
  baseline_hazard <- c(exp(phi))

  # Normalize the baseline hazard plot
  area <- 0
  for(i in 1:L){
    area <- area  + (time_axis[i+1] - time_axis[i]) * baseline_hazard[i]
  }
  baseline_hazard <- baseline_hazard / area

  # Return the baseline hazard step-function
  return(baseline_hazard)
}

