#' @title
#' Posterior Frailty Estimates
#'
#' @description
#' Function for computing the posterior frailty estimates of the time-dependent shared frailty Cox model.
#' 
#' Recalling the structure of the frailty \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k} 
#' with \eqn{k=1,\dots,L} and \eqn{j=1,\dots,N} as being composed by the sum
#' of two independent gamma distributions:
#' - \eqn{\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j}
#' - \eqn{\epsilon_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k}
#' 
#' The posterior frailty estimate is \eqn{\hat{Z}_{jk} = \hat{\alpha}_{j}/\hat{\alpha}_{max} + \hat{\epsilon}_{jk}/\hat{\epsilon}_{max}}. 
#' This function allows to get either the entire posterior frailty estimate \eqn{\hat{Z}_{jk}}
#' or its time-independent \eqn{\frac{\hat{\alpha}_{j}}{\hat{\alpha}_{\text{max}}}} or 
#' time-dependent \eqn{\frac{\hat{\epsilon}_{jk}}{\hat{\epsilon}_{\text{max}}}} components.
#' The user can control which components to display using the flag_eps and flag_alpha parameters. 
#' Only one of these flags can be set to TRUE at a time.
#'
#' @param object S3 object of class 'AdPaik' returned by the main model output, that contains all the information for the computation
#' of the frailty standard deviation.
#' @param flag_eps Logical flag indicating whether to extract only the time-dependent posterior frailty estimates. Default is FALSE.
#' @param flag_alpha Logical flag indicating whether to extract only the time-independent posterior frailty estimates. Default is FALSE.
#'
#' @return Vector or matrix of posterior frailty estimates, depending on the flag_eps and flag_alpha values. 
#' Specifically:
#'  - It is a vector of length equal to the N containing posterior frailty estimates for \eqn{\alpha_j, \forall j}.
#'    In this case the flag_eps must be FALSE and the flag_alpha must be TRUE.
#'  - Matrix of dimension (N, L) containing posterior frailty estimates for \eqn{\epsilon_{jk}, \forall j,k}.
#'    In this case the flag_eps must be TRUE and the flag_alpha must be FALSE.
#'  - Matrix of dimension (N, L) containing posterior frailty estimates for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. 
#'    In this case the flag_eps must be FALSE and the flag_alpha must be FALSE.
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
#' post_frailty_est(result)
#' }
#' 
post_frailty_est <- function(object, flag_eps = FALSE, flag_alpha = FALSE){
  
  if(flag_eps & flag_alpha)
    stop("Both flags cannot be TRUE at the same time.")
  else if(flag_eps & !flag_alpha){
    value = object$PosteriorFrailtyEstimates$eps
  } 
  else if(!flag_eps & flag_alpha){
    value = object$PosteriorFrailtyEstimates$alpha
  }
  else{
    value = object$PosteriorFrailtyEstimates$Z
  }
  return(value)
}



#_______________________________________________________________________________________________________________________


#' @title
#' Posterior Frailty Variances
#'
#' @description
#' Function for computing the posterior frailty variances of the time-dependent shared frailty Cox model.
#' 
#' Recalling the structure of the frailty \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k} 
#' with \eqn{k=1,\dots,L} and \eqn{j=1,\dots,N} as being composed by the sum
#' of two independent gamma distributions:
#' - \eqn{\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j}
#' - \eqn{\epsilon_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k}
#' 
#' The posterior frailty variance is \eqn{var(\hat{Z}_{jk}) = var(\hat{\alpha}_{j}/\hat{\alpha}_{max}) + var(\hat{\epsilon}_{jk}/\hat{\epsilon}_{max}}). 
#' This function allows to get either the entire posterior frailty variance \eqn{var(\hat{Z}_{jk})}
#' or its time-independent \eqn{var(\frac{\hat{\alpha}_{j}}{\hat{\alpha}_{\text{max}}})} or 
#' time-dependent \eqn{var(\frac{\hat{\epsilon}_{jk}}{\hat{\epsilon}_{\text{max}}})} components.
#' The user can control which components to display using the flag_eps and flag_alpha parameters. 
#' Only one of these flags can be set to TRUE at a time.
#'
#' @param object S3 object of class 'AdPaik' returned by the main model output, that contains all the information for the computation
#' of the frailty standard deviation.
#' @param flag_eps Logical flag indicating whether to extract only the time-dependent posterior frailty estimates. Default is FALSE.
#' @param flag_alpha Logical flag indicating whether to extract only the time-independent posterior frailty estimates. Default is FALSE.
#'
#' @return Vector or matrix of posterior frailty variances, depending on the flag_eps and flag_alpha values. 
#' Specifically:
#'  - It is a vector of length equal to the N containing posterior frailty variances for \eqn{\alpha_j, \forall j}.
#'    In this case the flag_eps must be FALSE and the flag_alpha must be TRUE.
#'  - Matrix of dimension (N, L) containing posterior frailty variances for \eqn{\epsilon_{jk}, \forall j,k}.
#'    In this case the flag_eps must be TRUE and the flag_alpha must be FALSE.
#'  - Matrix of dimension (N, L) containing posterior frailty variances for \eqn{Z_{jk} \forall j,k}. 
#'    In this case the flag_eps must be FALSE and the flag_alpha must be FALSE.
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
#' post_frailty_var(result)
#' }
#' 
#' @export
post_frailty_var <- function(object, flag_eps = FALSE, flag_alpha = FALSE){
  
  if(flag_eps & flag_alpha)
    stop("Both flags cannot be TRUE at the same time.")
  else if(flag_eps & !flag_alpha){
    value = object$PosteriorFrailtyVariance$epsVar
  } 
  else if(!flag_eps & flag_alpha){
    value = object$PosteriorFrailtyVariance$alphaVar
  }
  else{
    value = object$PosteriorFrailtyVariance$ZVar
  }
  return(value)
}



#_______________________________________________________________________________________________________________________


#' @title
#' Posterior Frailty Estimates and Variances for the 'Adapted Paik et Al.'s Model'
#'
#' @param optimal_params Optimal parameters estimated by maximizing the log-likelihood function, through the constraint
#' multi-dimensional optmization method.
#' @param dataset Dataset containing all the covariates/regressors.
#' @param time_to_event Time-instant, in the follow-up, in which an individual faces the event or fails.
#' If an individual does not face the event in the follow-up, then the time-instant must assume a default value.
#' @param centre Individual group/cluster membership.
#' @param time_axis Temporal domain.
#'
#' @return S3 object of class 'PF.AdPaik' composed of two elements of different class:
#' - PosteriorFrailtyEst: S3 object of class 'PFE.AdPaik'.
#' - PosteriorFrailtyVar: S3 object of class 'PFV.AdPaik'.
#'
#' @keywords internal

post_frailty_internal <- function(optimal_params, dataset, time_to_event, centre, time_axis){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  R <- n_regressors <- dim(dataset)[2]
  L <- n_intervals <- length(time_axis) - 1
  n_params <- length(optimal_params)
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Extract parameters from the vector
  phi <- optimal_params[1:L]
  betar <- optimal_params[(L + 1):(L + R)]
  mu1 <- optimal_params[(L + R + 1)]
  nu <- optimal_params[(L + R + 2)]
  gammak <- optimal_params[(L + R + 3) : n_params]
  
  # For idenfiability purpose
  mu2 <- 1 - mu1
  
  # Compute parameters for the gamma densities
  par1_alpha <- mu1/nu
  par2_alpha <- 1/nu
  par1_eps <- mu2/gammak
  par2_eps <- 1/gammak
  
  # Compute other dataset related variables
  extraction_variables <- extract_event_data(dataset, time_to_event, centre, time_axis, phi, betar)
  N_ik <- extraction_variables$N_ik
  N_i <- extraction_variables$N_i
  e_ijk <- extraction_variables$e_ijk
  Y_risk <- extraction_variables$Y_risk
  cum_hazard_group <- extraction_variables$cum_hazard_group
  sum_cum_hazard_group <- extraction_variables$sum_cum_hazard_group
  
  # Compute both frailty terms (alpha, eps) and their variance
  eps_frailty <- matrix(rep(0, n_intervals * n_centres), n_centres, n_intervals)
  eps_frailty_var <- eps_frailty
  for(i in 1:n_centres){
    for (k in 1:n_intervals) {
      eps_frailty[i,k] <- (par1_eps[k] + N_ik[i,k])/(par2_eps[k] + cum_hazard_group[i,k])
      eps_frailty_var[i,k] <- eps_frailty[i,k] / (par2_eps[k] + cum_hazard_group[i,k])
    }
  }
  
  alpha_frailty <- alpha_frailty_var <- rep(0, n_centres)
  for(i in 1:n_centres){
    alpha_frailty[i] <- (par1_alpha + N_i[i])/(par2_alpha + sum_cum_hazard_group[i])
    alpha_frailty_var[i] <- alpha_frailty[i] / (par2_alpha + sum_cum_hazard_group[i])
  }
  
  # Divide each term to have frailty mean equal to 1
  eps_max <- max(eps_frailty)
  alpha_max <- max(alpha_frailty)
  
  eps_frailty <- eps_frailty/eps_max
  alpha_frailty <- alpha_frailty/alpha_max
  
  # Sum both frailty terms to have the full estimated frailty
  Z_frailty <- eps_frailty + alpha_frailty
  
  post_frailty_est <- list("alpha" = alpha_frailty,
                           "eps" = eps_frailty,
                           "Z" = Z_frailty)
  class(post_frailty_est) <- "PFE.AdPaik"
  
  # Work on the variance to follow the previous relation:
  alpha_frailty_var <- alpha_frailty_var/(alpha_max^2)
  eps_frailty_var <- eps_frailty_var/(eps_max^2)
  Z_frailty_var <- eps_frailty_var + alpha_frailty_var
  
  post_frailty_var <- list("alphaVar" = alpha_frailty_var,
                           "epsVar" = eps_frailty_var,
                           "ZVar" = Z_frailty_var)
  class(post_frailty_var) <- "PFV.AdPaik"
  
  return_list <- list("PostFrailtyEst" = post_frailty_est,
                      "PostFrailtyVar" = post_frailty_var)
  class(return_list) <- "PF.AdPaik"
  
  return (return_list)
}

#-------------------------------------------------------------------------------
#' @title
#' Extracting Variables for Posterior Frailty Estimates Computation
#'
#' @description
#' Function for extracting from the dataset quantities necessary to the evaluation of the posterior
#' frailty estimates.
#'
#' @param dataset Dataset containing the covariates/regressors. Their numerosity is indicated with R.
#' @param time_to_event Time-instant in the follow-up in which an individual fails or faces the event.
#' If an individual does not face the event, the time-instant assumes a default value.
#' @param centre Categorical vector indicating the group/cluster membership. The number of distinct group is indicated with N.
#' @param time_axis Numerical vector of the temporal domain. Its length is (L+1), where L indicates the number of intervals of the time-domain.
#' @param phi Numerical vector of length L, of estimated baseline log-hazard.
#' @param betar Numerical vector of length R, of estimated regressors.
#'
#' @return S3 object of class 'EventData', composed of six elements. See details.
#'
#' @details
#' The S3 class obejct 'EventData' contains the variables necessary for the estimate of the posterior frailty and that can be extracted or
#' computed starting from the dataset.
#' - N_ik: matrix of dimension (N, L), containing the number of event in each interval k and group i.
#' - N_i: numerical vector of length L, with the number of event in each group i. It can be computed as: \eqn{N_i = \sum_{k=1}^L N_{ik}}.
#' - e_ijk: matrix of dimension (n_individuals, L) with the evaluation of the temporal integral, for each individual j, group i and interval k.
#' - Y_risk: binary matrix of dimension (n_individuals, L) reporting for each individual, in each interval, his/her risk of facing the event.
#' For an individual, the risk is equal to 1 in an interval k if, in that interval, he/she has not faced the event yet; otherwise, it is equal to 0.
#' - cum_hazard_group: matrix of dimension (N, L), where each element in position (i,k) indicates the computed cumulative hazard for all individuals
#' belonging to group i and at interval k.
#' - sum_cum_hazard_group: numerical vector of length N, giving the sum of the computed cumulative hazard for all intervals k and for all individuals
#' belonging to group i.
#' It can be computed from the previous element, summing with respect to the interval k.
#' 
#' @keywords internal

extract_event_data <- function(dataset, time_to_event, centre, time_axis, phi, betar){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Compute the number of dropout students in each interval
  N_ik <- matrix(rep(0, n_intervals * n_centres), n_centres, n_intervals)
  for(i in 1:n_centres){
    indexes_centre <- which(centre == centre_codes[i])
    n_individuals_centre <- length(indexes_centre)
    counter <- rep(0, n_intervals)
    for(j in 1:n_individuals_centre){
      index_j <- indexes_centre[j]
      for(k in 1:n_intervals){
        if(k != n_intervals){
          if((time_to_event[index_j] < time_axis[k+1]) & (time_to_event[index_j] >= time_axis[k]))
            counter[k] <- counter[k] + 1
        }
        else{
          if((time_to_event[index_j] <= time_axis[k+1]) & (time_to_event[index_j] >= time_axis[k]))
            counter[k] <- counter[k] + 1
        }
      }
    }
    N_ik[i,] <- counter
  }
  
  # Compute the number of dropout student in each faculty
  N_i <- rep(0, n_centres)
  for(i in 1: n_centres){
    indexes_centre <- which(centre == centre_codes[i])
    counter <- 0
    for(ii in 1:length(indexes_centre)){
      index <- indexes_centre[ii]
      if(time_to_event[index] < max(time_to_event))
        counter <- counter + 1
    }
    N_i[i] <- counter
  }
  
  # Define the variable e_ijk
  e_ijk <- matrix(rep(0, n_individuals * n_intervals), n_individuals, n_intervals)
  for(j in 1:n_individuals){
    for(k in 1:n_intervals){
      e_ijk[j,k] <- time_int_eval(time_to_event[j], k, time_axis)
    }
  }
  
  # Define the at risk indicator for each individual and interval
  Y_risk <- matrix(rep(0, n_individuals * n_intervals), n_individuals, n_intervals)
  for(j in 1:n_individuals){
    for(k in 1:n_intervals){
      if((time_to_event[j] < max(time_to_event)) & (time_to_event[j] > time_axis[k+1]))
        Y_risk[j,k] <- 1
      else if(time_to_event[j] == max(time_to_event))
        Y_risk[j,k] <- 1
    }
  }
  
  # Define the cumulative hazard for all individuals
  cum_hazard_jk <- matrix(rep(0, n_individuals * n_intervals), n_individuals, n_intervals)
  for(j in 1:n_individuals){
    dataset_betar <- as.numeric(dataset[j,]) %*% betar
    for(k in 1:n_intervals){
      cum_hazard_jk[j,k] <- e_ijk[j,k] * Y_risk[j,k] * exp(phi[k] + dataset_betar)
    }
  }
  
  # Define the cumulative hazard for each group
  cum_hazard_group <- matrix(rep(0, n_centres * n_intervals), n_centres, n_intervals)
  for(i in 1:n_centres){
    indexes_centre <- which(centre == centre_codes[i])
    n_individuals_group <- length(indexes_centre)
    for(k in 1:n_intervals)
      cum_hazard_group[i,k] <- sum(cum_hazard_jk[indexes_centre,k])
  }
  
  # Define the cumulative hazard for each faculty, time-independent
  sum_cum_hazard_group <- rep(0, n_centres)
  for(i in 1:n_centres){
    sum_cum_hazard_group[i] <- sum(cum_hazard_group[i,])
  }
  
  return_list <- list("N_ik" = N_ik,
                      "N_i" = N_i,
                      "e_ijk" = e_ijk,
                      "Y_risk" = Y_risk,
                      "cum_hazard_group" = cum_hazard_group,
                      "sum_cum_hazard_group" = sum_cum_hazard_group)
  class(return_list) <- "EventData"
  
  return (return_list)
}
#-------------------------------------------------------------------------------
#' @title
#' Confidence Interval for Posterior Frailty Estimates
#'
#' @description
#' Function for computing the confidence interval for each posterior frailty estimates \eqn{\hat{Z}_{jk}}.
#'
#' @param post_frailty_est Posterior frailty estimates list.
#' @param post_frailty_est_var Posterior frailty variance list.
#' @param n_centres Number of clusters/centres.
#' @param n_intervals Number of intervals of the time-domain. it is equal to the length of the tima_axis minus one.
#' @param level A numeric value representing the confidence level.
#'
#' @return S3 object of class 'PFCI.AdPaik' composed of two matrices of dimension (number groups, number of intervals):
#' - PostFrailtyCI_left: left confidence interval for each posterior frailty estimates
#' - PostFrailtyCI_right: right confidence interval for each each posterior frailty estimates
#' 
#' @keywords internal

post_frailty_CI_internal <- function(post_frailty_est, post_frailty_est_var, n_centres, n_intervals, level){
  # Check structure correctness
  check.structure_post_frailty_est(post_frailty_est, n_intervals, n_centres)
  
  # Check value of posterior frailty estimates
  check.value_post_frailty(post_frailty_est, n_centres, n_intervals)
  check.value_post_frailty(post_frailty_est_var, n_centres, n_intervals)
  
  post_frailty_CI_left <- matrix(rep(0, n_intervals * n_centres), n_centres, n_intervals)
  post_frailty_CI_right <-  matrix(rep(0, n_intervals * n_centres), n_centres, n_intervals)
  
  alpha <- 1 - level
  z_critical <- stats::qnorm(1 - alpha / 2)  # Two-tailed
  
  for(i in 1:n_centres){
    for(k in 1:n_intervals){
      sd_est <- sqrt(post_frailty_est_var$ZVar[i,k])
      post_frailty_CI_left[i,k] <- post_frailty_est$Z[i,k] - z_critical * sd_est
      post_frailty_CI_right[i,k] <- post_frailty_est$Z[i,k] + z_critical * sd_est
    }
  }
  post_frailty_CI <- list("PostFrailtyCI_left" = post_frailty_CI_left,
                          "PostFrailtyCI_right" = post_frailty_CI_right)
  class(post_frailty_CI) <- "PFCI.AdPaik"
  
  return (post_frailty_CI)
}



#_______________________________________________________________________________________________________________________


#' @title
#' Posterior Frailty Confidence Intervals
#'
#' @description
#' Function for computing the posterior frailty confidence intervals of the time-dependent shared frailty Cox model.
#' 
#' Recalling the structure of the frailty \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k} 
#' with \eqn{k=1,\dots,L} and \eqn{j=1,\dots,N} as being composed by the sum
#' of two independent gamma distributions:
#' - \eqn{\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j}
#' - \eqn{\epsilon_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k}
#' 
#' The posterior frailty estimate is \eqn{\hat{Z}_{jk} = \hat{\alpha}_{j}/\hat{\alpha}_{max} + \hat{\epsilon}_{jk}/\hat{\epsilon}_{max}}. 
#' This function allows to get the confidence intervals of either the entire posterior frailty estimates \eqn{\hat{Z}_{jk}}
#' or its time-independent \eqn{\frac{\hat{\alpha}_{j}}{\hat{\alpha}_{\text{max}}}} or 
#' time-dependent \eqn{\frac{\hat{\epsilon}_{jk}}{\hat{\epsilon}_{\text{max}}}} components.
#' The user can control which components to display using the flag_eps and flag_alpha parameters. 
#' Only one of these flags can be set to TRUE at a time.
#'
#' @param object S3 object of class 'AdPaik' returned by the main model output, that contains all the information for the computation
#' of the frailty standard deviation.
#' @param level A numeric value representing the confidence level for the posterior frailty confidence intervals.
#' Default is 0.95 for 95% confidence.
#' @param flag_eps Logical flag indicating whether to extract only the time-dependent posterior frailty estimates. Default is FALSE.
#' @param flag_alpha Logical flag indicating whether to extract only the time-independent posterior frailty estimates. Default is FALSE.
#'
#' @return A list for posterior frailty confidence intervals, depending on the flag_eps and flag_alpha values. 
#' Specifically:
#'  - A list of length equal to the N containing posterior frailty confidence intervals for \eqn{\alpha_j, \forall j}.
#'    In this case the flag_eps must be FALSE and the flag_alpha must be TRUE.
#'  - A list of length equal to the NxL containing posterior frailty confidence intervals for \eqn{\epsilon_{jk}, \forall j,k}.
#'    In this case the flag_eps must be TRUE and the flag_alpha must be FALSE.
#'  - A list of length equal to the NxL containing posterior frailty confidence intervals for \eqn{Z_{jk} \forall j,k}. 
#'    In this case the flag_eps must be FALSE and the flag_alpha must be FALSE.
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
#' post_frailty_confint(result)
#' }
#' 
#' @export
post_frailty_confint <- function(object, level = 0.95, flag_eps = FALSE, flag_alpha = FALSE) {
  
  alpha <- 1 - level
  z_critical <- stats::qnorm(1 - alpha / 2)  # Two-tailed
  
  n_centres = object$NClusters
  n_intervals = object$NIntervals
  
  # Format confidence levels
  a <- (1 - level) / 2
  pct <- format_perc(c(a, 1 - a), 3)
  
  if (flag_alpha) {
    # Handle alpha confidence intervals
    alpha_est <- object$PosteriorFrailtyEstimates$alpha
    alpha_sd <- sqrt(object$PosteriorFrailtyVariance$alphaVar)
    
    alpha_CI_left <- alpha_est - z_critical * alpha_sd
    alpha_CI_right <- alpha_est + z_critical * alpha_sd
    
    # Create named vector with confidence intervals
    alpha_ci <- cbind(alpha_CI_left, alpha_CI_right)
    colnames(alpha_ci) <- pct
    rownames(alpha_ci) <- object$ClusterCodes
    
    return(alpha_ci)
    
  } else {
    # Initialize matrices
    post_frailty_CI_left <- matrix(0, nrow = n_centres, ncol = n_intervals)
    post_frailty_CI_right <- matrix(0, nrow = n_centres, ncol = n_intervals)
    
    for (i in 1:n_centres) {
      for (k in 1:n_intervals) {
        if(flag_eps == FALSE){
          sd_est <- sqrt(object$PosteriorFrailtyVariance$ZVar[i, k])
          post_frailty_CI_left[i, k] <- object$PosteriorFrailtyEstimates$Z[i, k] - z_critical * sd_est
          post_frailty_CI_right[i, k] <- object$PosteriorFrailtyEstimates$Z[i, k] + z_critical * sd_est
        }
        else if(flag_eps == TRUE){
          sd_est <- sqrt(object$PosteriorFrailtyVariance$epsVar[i, k])
          post_frailty_CI_left[i, k] <- object$PosteriorFrailtyEstimates$eps[i, k] - z_critical * sd_est
          post_frailty_CI_right[i, k] <- object$PosteriorFrailtyEstimates$eps[i, k] + z_critical * sd_est
        }
      }
    }
    
    # Combine into a single matrix
    ci_matrix <- matrix(NA_real_, nrow = n_centres * n_intervals, ncol = 2)
    ci_matrix[, 1] <- as.vector(post_frailty_CI_left)
    ci_matrix[, 2] <- as.vector(post_frailty_CI_right)
    
    # Assign row and column names
    rownames(ci_matrix) <- paste0(rep(object$ClusterCodes, times = n_intervals), 
                                  "_Interval_", rep(1:n_intervals, each = n_centres))
    colnames(ci_matrix) <- pct
    
    return(ci_matrix)
    
  }

}




