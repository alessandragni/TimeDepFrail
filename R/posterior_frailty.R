#' @title
#' Posterior frailty estimates and variances for the 'Adapted Paik et al.'s Model'
#'
#' @description
#' Function for computing the posterior frailty estimates and variances of the time-dependent shared frailty Cox model.
#' Recalling the structure of the frailty \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k} as being composed by the sum
#' of two independent gamma distributions:
#' - \eqn{\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j}
#' - \eqn{\epsilon_{jk} \sin gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k}
#' the posterior distribution of both terms is still a gamma with different mean and variance and the
#' posterior frailty estimate corresponds to the 'empirical Bayes estimate', that is the previous mentioned posterior mean.
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

post_frailty.AdPaik <- function(optimal_params, dataset, time_to_event, centre, time_axis){
  
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
#' Extracting variables for Posterior Frailty Estimates computation
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
      if(time_to_event[index] < 6.1)
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
      if((time_to_event[j] < 6.1) & (time_to_event[j] > time_axis[k+1]))
        Y_risk[j,k] <- 1
      else if(time_to_event[j] == 6.1)
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
#' Confidence interval for posterior frailty estimates
#'
#' @description
#' Function for computing the confidence interval for each posterior frailty estimates \eqn{\hat{Z}_{jk}}.
#'
#' @param post_frailty_est Posterior frailty estimates list.
#' @param post_frailty_est_var Posterior frailty variance list.
#' @param n_centres Number of clusters/centres
#' @param n_intervals Number of intervals of the time-domain. it is equal to the length of the tima_axis minus one.
#'
#' @return S3 object of class 'PFCI.AdPaik' composed of two matrices of dimension (number groups, number of intervals):
#' - PostFrailtyCI_left: left confidence interval for each posterior frailty estimates
#' - PostFrailtyCI_right: right confidence interval for each each posterior frailty estimates

post_frailty_CI.AdPaik <- function(post_frailty_est, post_frailty_est_var, n_centres, n_intervals){
  # Check structure correctness
  check.structure_post_frailty_est(post_frailty_est, n_intervals, n_centres)
  
  # Check value of posterior frailty estimates
  check.value_post_frailty(post_frailty_est, n_centres, n_intervals)
  check.value_post_frailty(post_frailty_est_var, n_centres, n_intervals)
  
  post_frailty_CI_left <- matrix(rep(0, n_intervals * n_centres), n_centres, n_intervals)
  post_frailty_CI_right <-  matrix(rep(0, n_intervals * n_centres), n_centres, n_intervals)
  for(i in 1:n_centres){
    for(k in 1:n_intervals){
      sd_est <- sqrt(post_frailty_est_var$ZVar[i,k])
      post_frailty_CI_left[i,k] <- post_frailty_est$Z[i,k] - 1.96 * sd_est
      post_frailty_CI_right[i,k] <- post_frailty_est$Z[i,k] + 1.96 * sd_est
    }
  }
  post_frailty_CI <- list("PostFrailtyCI_left" = post_frailty_CI_left,
                          "PostFrailtyCI_right" = post_frailty_CI_right)
  class(post_frailty_CI) <- "PFCI.AdPaik"
  
  return (post_frailty_CI)
}
