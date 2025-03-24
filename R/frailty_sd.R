#-------------------------------------------------------------------------------

#' @title
#' Frailty Standard Deviation and Variance for the 'Adapted Paik et Al.'s Model'
#'
#' @description
#' The function computes both the standard deviation and the variance of the time-dependent
#' frailty of the 'Adapted Paik et al.'s Model'.
#'
#' Recalling the frailty structure \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}} as being composed by a constant group-dependent term
#' (\eqn{\alpha_j}) and a time and group dependent term (\eqn{\epsilon_{jk}}), the frailty variance (and standard deviation)
#' can be computed in two different way:
#' - Considering only the time-dependent spread of the clusters/groups/centre: \eqn{var(Z_{jk}) = \mu_2 * \gamma_k}.
#' In this case, the flag_full should be FALSE and flag_variance should be TRUE.
#'
#' - Considering both the time-dependent and constant spread of the clusters: \eqn{var(Z_{jk}) = \mu_1 * \nu + \mu_2 * \gamma_k}.
#' The new added term only moves upward the other case and the flag_full should be TRUE and flag_variance should be TRUE.
#'
#' @param object S3 object of class 'AdPaik' returned by the main model output, that contains all the information for the computation
#' of the frailty standard deviation.
#' @param flag_full A boolean flag indicating whether to get the full standard deviation (`TRUE`) or only the time-dependent component (`FALSE`). Default to `TRUE`.
#' @param flag_variance A boolean flag indicating whether to get the frailty variance (`TRUE`) or the frailty standard deviation (`FALSE`). Default to `FALSE`.
#'
#' @return Numerical vector of length equal to the number of intervals of the time-domain, 
#' with the value of the frailty standard deviation or variance (either full or only the time-dependent component).
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
#' frailty_sd(result)
#' }
frailty_sd <- function(object, flag_full = TRUE, flag_variance = FALSE) {
  
  # Extract information from the object
  optimal_params <- object$OptimalParameters
  time_axis <- object$TimeDomain
  n_regressors <- object$NRegressors
  categories_range_min <- object$ParametersRange$ParametersRangeMin
  categories_range_max <- object$ParametersRange$ParametersRangeMax
  
  frailty_Sd_internal = frailty_Sd_internal(optimal_params, time_axis, n_regressors, 
                                            categories_range_min, categories_range_max, flag_full)
  
  if (flag_variance)
    return(frailty_Sd_internal$FrailtyVariance)
  else
    return(frailty_Sd_internal$FrailtyStandardDeviation)
  
}

#-------------------------------------------------------------------------------

#' @title
#' Internal Function for Frailty Standard Deviation for the 'Adapted Paik et Al.'s Model'
#'
#' @description
#' The function computes both the standard deviation and the variance of the time-dependent
#' frailty of the 'Adapted Paik et al.'s Model'.
#'
#' @param optimal_params Optimal parameter vector, estimated through multi-dimensional optimization of the log-likelihood function.
#' @param time_axis Partition of the temporal domain.
#' @param n_regressors Number of regressors of the dataset. This value must be provided to be able to correctly
#' compute the number of parameters of the model and, therefore, to check the dimension of the parameter vector.
#' @param categories_range_min Vector of minimum value (range) assumed by the parameters category.
#' @param categories_range_max Vector of maximum value (range) assumed by the parameters category.
#' @param flag_full Do we want to compute the full frailty standard deviation (second case)? If so, the flag must be TRUE,
#' otherwise (first case), FALSE.
#'
#' @return S3 class object 'FrailtyDispersion' containing both two numerical vectors of length equal to the number of intervals of the time-domain:
#' - FrailtyVariance
#' - FrailtyStandardDeviation
#'
#' @keywords internal
frailty_Sd_internal <- function(optimal_params, time_axis, n_regressors,
                                categories_range_min, categories_range_max,
                                flag_full) {
  
  # Extract information from input variables
  L <- n_intervals <- length(time_axis) - 1
  R <- n_regressors
  
  # Extract parameters from optimal vector
  mu1 <- optimal_params[L + R + 1]
  nu <- optimal_params[(L + 2 + R)]
  gammak <- optimal_params[(L + 3 + R):(2 * L + R + 2)]
  
  # For idenfiability purpose
  mu2 <- 1 - mu1
  
  # Compute the frailty variance and standard deviation
  variance <- sd <- rep(0, L)
  variance_k <- 0
  for (k in 1:L) {
    if (flag_full)
      variance_k <- mu1 * nu + mu2 * gammak[k]
    else
      variance_k <- mu2 * gammak[k]
    
    if (variance_k < 0) {
      msg <- paste('Negative frailty variance in position ', k, '.')
      stop(msg)
    }
    
    variance[k] <- variance_k
    sd[k] <- sqrt(variance[k])
  }
  
  return_list <- list("FrailtyVariance" = variance,
                      "FrailtyStandardDeviation" = sd)
  class(return_list) <- 'FrailtyDispersion'
  
  return(return_list)
}


