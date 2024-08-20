#' @title
#' Estimated frailty mean for the 'Centre-Specific Frailty Model with Power Parameter'
#'
#' @description
#' The function computes the model estimated frailty mean, giving the estimated parameters.
#'
#' @param optimal_params Numerical vector of optimal estimated parameters, of length equal to the number of
#' model parameters.
#' @param time_axis Numerical vector of the time-domain.
#' @param n_regressors Number of regressors of the model
#'
#' @return Numerical vector of length equal to the number of intervals, containing the estimated frailty 
#' mean of each interval.

frailty_mean.PowPar <- function(optimal_params, time_axis, n_regressors){
  # Extract information from input variables
  L <- length(time_axis) - 1
  R <- n_regressors

  # Extract parameters from the optimal vector
  gammak <- rep(0,L)
  gammak[2:L] <- optimal_params[(L+1+R):(2*L+R-1)]
  gammak[1] <- 1
  sigma <- optimal_params[2*L+R]

  # Compute estimated frailty mean
  frailty_mean <- rep(0,L)
  for(k in 1:L){
    frailty_mean[k] <- exp(((sigma * gammak[k])^2)/2)
  }

  return (frailty_mean)
}

