#' @title
#' Baseline hazard step-function
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

bas_hazard <- function(optimal_params, time_axis){
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
