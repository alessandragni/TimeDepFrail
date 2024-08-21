# #' @title
# #' Frailty standard deviation for the 'Adapted Paik et al.'s Model'
# #'
# #' @description
# #' The function computes both the standard deviation and the variance of the time-dependent
# #' frailty of the 'Adapted Paik et al.'s Model'.
# #'
# #' Recalling the frailty structure \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}} as being composed by a constant group-dependent term
# #' (\eqn{\alpha_j}) and a time and group dependent term (\eqn{\epsilon_{jk}}), the frailty standard deviation (and variance)
# #' can be computed in two different way:
# #' - Considering only the time-dependent spread of the clusters/groups/centre: \eqn{sd(Z_{jk}) = \mu_2 * \gamma_k}.
# #' In this case, the flag_fullsd should be FALSE.
# #'
# #' - Considering both the time-dependent and constant spread of the clusters: \eqn{sd(Z_{jk}) = \mu_1 * \nu \mu_2 * \gamma_k}.
# #' The new added term only moves upward the other case and the flag_fullsd should be TRUE.
# #'
# #' The final case only depends on what we want to observe.
# #'
# #' @param optimal_params Optimal parameter vector, estimated through multi-dimensional optimization of the log-likelihood function.
# #' @param time_axis Partition of the temporal domain.
# #' @param n_regressors Number of regressors of the dataset. This value must be provided to be able to correctly
# #' compute the number of parameters of the model and, therefore, to check the dimension of the parameter vector.
# #' @param categories_range_min Vector of minimum value (range) assumed by the parameters category.
# #' @param categories_range_max Vector of maximum value (range) assumed by the parameters category.
# #' @param flag_fullsd Do we want to compute the full frailty standard deviation (second case)? If so, the flag must be TRUE,
# #' otherwise (first case), FALSE.
# #'
# #' @return S3 class object 'FrailtyDispersion' containing both two numerical vectors of length equal to the numbero of intervals of the time-domain:
# #' - FrailtyVariance
# #' - FrailtyStandardDevation
# 
# frailty_sd.AdPaik <- function (optimal_params, time_axis, n_regressors,
#                         categories_range_min, categories_range_max,
#                         flag_fullsd = TRUE){
# 
#   # Extract information from input variables
#   L <- n_intervals <- length(time_axis) - 1
#   R <- n_regressors
#   n_params <- length(optimal_params)
# 
#   # Define vector of categories for Adapted Paik et al.'s Model
#   params_categories <- c(n_intervals, n_regressors, 1, 1, n_intervals)
#   n_categories <- length(params_categories)
# 
#   # Check correctness of input categories
#   check.categories_params(n_categories, categories_range_min, categories_range_max)
# 
#   # Check correctness of input optimal parameter vector
#   if(n_params != (2*n_intervals + n_regressors + 2))
#     stop("Provided 'optimal_params' vector of length different from theoretical one for current model.")
# 
#   # Generate extended vector of parameters ranges
#   params_range_min <- params_range_max <- c()
#   for(c in 1: n_categories){
#     n_params_in_c <- params_categories[c]
#     params_range_min <- c(params_range_min, rep(categories_range_min[c], n_params_in_c))
#     params_range_max <- c(params_range_max, rep(categories_range_max[c], n_params_in_c))
#   }
# 
#   # Controll optimal_parameters are contained in the min and max range
#   check.range_params(optimal_params, params_range_min, params_range_max)
# 
#   # Extract parameters from optimal vector
#   mu1 <- optimal_params[L+R+1]
#   nu <- optimal_params[(L+2+R)]
#   gammak <- optimal_params[(L+3+R):(2*L+R+2)]
# 
#   # For idenfiability purpose
#   mu2 <- 1 - mu1
# 
#   # Compute the frailty variance and standard deviation
#   variance <- sd <- rep(0, L)
#   variance_k <- 0
#   for (k in 1:L){
#     if(flag_fullsd)
#       variance_k <- mu1 * nu + mu2 * gammak[k]
#     else
#       variance_k <- mu2 * gammak[k]
#     
#     if(variance_k < 0){
#       msg <- paste('Negative frailty variance in position ', k, '.')
#       stop(msg)
#     }
# 
#     variance[k] <- variance_k
#     sd[k] <- sqrt(variance[k])
#   }
#   
#   return_list <- list("FrailtyVariance" = variance,
#                       "FrailtyStandardDeviation" = sd)
#   class(return_list) <- 'FrailtyDispersion'
# 
#   return (return_list)
# }
#-------------------------------------------------------------------------------
#' @title
#' Frailty standard deviation and Variance for the 'Adapted Paik et al.'s Model'
#'
#' @description
#' The function computes both the standard deviation and the variance of the time-dependent
#' frailty of the 'Adapted Paik et al.'s Model'.
#'
#' Recalling the frailty structure \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}} as being composed by a constant group-dependent term
#' (\eqn{\alpha_j}) and a time and group dependent term (\eqn{\epsilon_{jk}}), the frailty standard deviation (and variance)
#' can be computed in two different way:
#' - Considering only the time-dependent spread of the clusters/groups/centre: \eqn{sd(Z_{jk}) = \mu_2 * \gamma_k}.
#' In this case, the flag_fullsd should be FALSE.
#'
#' - Considering both the time-dependent and constant spread of the clusters: \eqn{sd(Z_{jk}) = \mu_1 * \nu \mu_2 * \gamma_k}.
#' The new added term only moves upward the other case and the flag_fullsd should be TRUE.
#'
#' The final case only depends on what we want to observe.
#'
#' @param result S3 object of class 'AdPaik' returned by the main model output, that contains all the information for the computation
#' of the frailty standard deviation.
#' @param flag_fullsd Logical value. Do we want to compute the full frailty standard deviation? If so, the flag must be TRUE,
#' otherwise, FALSE.
#'
#' @return S3 class object 'FrailtyDispersion' containing both two numerical vectors of length equal to the numbero of intervals of the time-domain:
#' - FrailtyVariance
#' - FrailtyStandardDevation
#'
#' @export
#'
#' @examples
#' # Define the variables needed for the model execution
#' formula <- time_to_event ~ GenderF + CFUP + cluster(group)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' eps <- 1e-10
#' categories_range_min <- c(-8, -2, eps, eps, eps)
#' categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)
#'
#' # Call the main model
#' result <- AdPaikModel(formula, data, time_axis,
#'                       categories_range_min, categories_range_max, TRUE)
#'
#' frailty_sd(result, TRUE)
#' frailty_sd(result, FALSE)

frailty_sd.AdPaik <- function (result, flag_fullsd = TRUE){

  # Check result structure
  check.result(result)

  # Extract information from input variables
  L <- n_intervals <- result$NIntervals
  R <- n_regressors <- result$NRegressors
  n_params <- result$NParameters
  optimal_params <- result$OptimalParameters

  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- result$ParametersCategories
  n_categories <- length(params_categories)

  # Generate extended vector of parameters ranges
  params_range_min <- result$ParametersRange$ParametersRangeMin
  params_range_max <- result$ParametersRange$ParametersRangeMax

  # Controll optimal_parameters are contained in the min and max range
  check.range_params(optimal_params, params_range_min, params_range_max)

  # Extract parameters from optimal vector
  mu1 <- optimal_params[L+R+1]
  nu <- optimal_params[(L+2+R)]
  gammak <- optimal_params[(L+3+R):(2*L+R+2)]

  # For idenfiability purpose
  mu2 <- 1 - mu1

  # Compute the frailty variance and standard deviation
  variance <- sd <- rep(0, L)
  variance_k <- 0
  for (k in 1:L){
    if(flag_fullsd)
      variance_k <- mu1 * nu + mu2 * gammak[k]
    else
      variance_k <- mu2 * gammak[k]

    if(variance_k < 0){
      msg <- paste("Negative frailty variance in position ", k, ".")
      stop(msg)
    }
    
    variance[k] <- variance_k
    sd[k] <- sqrt(variance[k])
  }
  
  return_list <- list("FrailtyVariance" = variance,
                      "FrailtyStandardDeviation" = sd)
  class(return_list) <- 'FrailtyDispersion'

  return (return_list)
}
#-------------------------------------------------------------------------------
frailty_sd.PowPar <- function (optimal_params, time_axis, n_regressors,
                               categories_range_min, categories_range_max){
  
  # Extract information from input variables
  L <- n_intervals <- length(time_axis) - 1
  R <- n_regressors
  n_params <- length(optimal_params)
  
  # Define vector of categories for Centre-Specific Frailty Model with Power Parameter
  params_categories <- c(n_intervals, n_regressors, n_intervals - 1, 1)
  n_categories <- length(params_categories)
  
  # Check correctness of input categories
  check.categories_params(n_categories, categories_range_min, categories_range_max)
  
  # Check correctness of input optimal parameter vector
  if(n_params != (2*n_intervals + n_regressors))
    stop("Provided 'optimal_params' vector of length different from theoretical one for current model.")
  
  # Generate extended vector of parameters ranges
  params_range_min <- params_range_max <- c()
  for(c in 1: n_categories){
    n_params_in_c <- params_categories[c]
    params_range_min <- c(params_range_min, rep(categories_range_min[c], n_params_in_c))
    params_range_max <- c(params_range_max, rep(categories_range_max[c], n_params_in_c))
  }
  
  # Controll optimal_parameters are contained in the min and max range
  check.range_params(optimal_params, params_range_min, params_range_max)
  
  # Extract parameters from the optimal vector
  gammak <- matrix(rep(0,L), nrow = L, ncol = 1)
  gammak[2:L,1] <- optimal_params[(L+1+R):(2*L+R-1)]                # frailty time-dependence parameter. The first element is assigned later
  gammak[1,1] <- 1
  sigma <- optimal_params[2*L+R]                                    # standard deviation of normal Yi  
  
  # Store the variance of the frailty
  variance <- sd <- rep(0, L)
  for (k in 1:L){
    variance[k] <- (sigma * gammak[k,1])^2
    sd[k] <- sqrt(variance[k])
  } #   exp(2*(gammak[k,1]*sigma)^2) - exp((gammak[k,1]*sigma)^2)      (sigma * gammak[k,1])^2
  
  return_list <- list("FrailtyVariance" = variance,
                      "FrailtyStandardDeviation" = sd)
  class(return_list) <- 'FrailtyDispersion'

  return (return_list)
}


#-------------------------------------------------------------------------------
frailty_sd.StocTimeDep <- function(optimal_params, time_domain, n_regressors,
                                   categories_range_min, categories_range_max){

  # Extract information from input variables
  L <- n_intervals <- length(time_axis) - 1
  R <- n_regressors
  n_params <- length(optimal_params)

  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, 1, 1, 1)
  n_categories <- length(params_categories)

  # Check correctness of input categories
  check.categories_params(n_categories, categories_range_min, categories_range_max)

  # Check correctness of input optimal parameter vector
  if(n_params != (n_intervals + n_regressors + 3))
    stop("Provided 'optimal_params' vector of length different from theoretical one for current model.")

  # Generate extended vector of parameters ranges
  params_range_min <- params_range_max <- c()
  for(c in 1: n_categories){
    n_params_in_c <- params_categories[c]
    params_range_min <- c(params_range_min, rep(categories_range_min[c], n_params_in_c))
    params_range_max <- c(params_range_max, rep(categories_range_max[c], n_params_in_c))
  }

  # Control optimal_parameters are contained in the min and max range
  check.range_params(optimal_params, params_range_min, params_range_max)

  # Extract information from input variables
  n_time_points <- length(time_domain)

  # Extract parameters
  lambda1 <- params[(n_params-2):(n_params-2)]
  lambda2 <- params[(n_params-1):(n_params-1)]
  angle_alpha <- params[(n_params):(n_params)]

  # Convert parameters
  sigma2c <- lambda1*(cos(angle_alpha))^2 + lambda2*(sin(angle_alpha))^2
  sigmacb <- (lambda1 - lambda2)*(sin(angle_alpha))*(cos(angle_alpha))
  sigma2b <- lambda1*(sin(angle_alpha))^2 + lambda2*(cos(angle_alpha))^2

  # Initialize the variance and standard deviation
  variance <- sd <- rep(0, n_time_points)
  for(t in 1:n_time_points){
    time <- time_domain[t]
    variance[t] <- sigma2c + sigma2b*(time^2) + 2*sigmacb*time
    sd[t] <- sqrt(variance[t])
  }

  return (list("FrailtyVariance" = variance,
               "FrailtyStandardDeviation" = sd))
}
#-------------------------------------------------------------------------------
#' @title Function for computing frailty standard deviation and Variance 
#' 
#' @description According to the model it has been called, 
#' this function inspects the class output and then call the proper summary function.
#'
#' @result Output of the frailty_sd call. The result has its own class, associated to the called model.

frailty_sd <- function(result, flag_fullsd = TRUE){
  if(class(result) == "AdPaik")
    frailty_sd.AdPaik(result, flag_fullsd = TRUE)
  else if(class(result) == "PowPar")
    frailty_sd.PowPar(result$OptimalParameters, result$TimeDomain, result$NRegressors,
                      result$CategoriesRangeMin, result$CategoriesRangeMax)
  else if(class(result) == "StocTimeDep")
    frailty_sd.StocTimeDep(result$OptimalParameters, result$TimeDomain, result$NRegressors,
                           result$CategoriesRangeMin, result$CategoriesRangeMax)
}

#-------------------------------------------------------------------------------


