#-------------------------------------------------------------------------------
#' @title
#' One-Dimensional Log-Likelihood Function to be Optimized
#'
#' @description
#' Model log-likelihood function to be optimized only with respect to a parameter. To correctly identify this parameter inside the model
#' and inside the vector of all parameter, it is necessary to provide also the position (index) of this parameter in the vector.
#'
#' This function is internally used by the main function @AdPaikModel to perform, as said, the one-dimensional optimization
#' through 'optimize'.
#' It cannot be used to evaluate the log-likelihood function at a vector of parameter and at the provided data. For this purpose, we
#' have to use another implemented function, called @ll_AdPaik_eval.
#'
#' @details
#' This function firstly divides the individuals according to their group/cluster membership, extracting group customized dataset and
#' other variables, and then compute the group log-likelihood function through the function @ll_AdPaik_centre_1D.
#' The produced group log-likelihood value is summed together the
#' other values into a unique result, that corresponds to the overall (and final) log-likelihood value.
#'
#' @param x Value of the parameter, with respect to which the log-likelihood function has to be optimized.
#' @param index Index of the parameter inside the parameter vector.
#' For instance, if we need to optimize the log-likelihood function with respect to the first regressor,
#' then @x will be generic but @index will be equal to (n_intervals + 1) because in the parameter vector
#' the first regressor appears after the baseline log-hazard group (n_intervals elements).
#' @param params Parameter vector.
#' @param dataset Matrix containing only the formula regressors, that is the regressors appearing in the
#' formula object provided by the user and eventually modified if they are categorical (nd therefore transformed into dummy variables).
#' @param centre Individual membership to the clusters.
#' @param time_axis Temporal domain.
#' @param dropout_matrix Binary matrix indicating in which interval of the time domain an individual failed. For an individual,
#' the sum of the row elements must be equal to 1 (if he/she failed) or 0 (if he/she does not failed).
#' It has dimension equal to (n_individuals, n_intervals).
#' @param e_matrix Matrix of dimension (n_individual, n_intervals), where each element contains the evaluation of the temporal
#' integral, performed through the function @param time_int_eval.
#'
#' @return Overall log-likelihood function
#' 
#' @keywords internal

ll_AdPaik_1D <- function(x, index, params, dataset, centre,
                         time_axis, dropout_matrix, e_matrix){
  
  # Extract information from the input variables
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- length(params)
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Define a variable that contains the log-likelihood result
  # Then update it by summing the partial result of each centre (loglik_centre)
  ll_overall <- 0
  
  # Loop over each centre and compute the log-likelihood associated to it
  for (i in 1:n_centres){
    # Extract individuals in a centre
    indexes_centre <- which(centre == centre_codes[i])
    
    # Check number of individuals in each centre
    if(length(indexes_centre) <= 1)
      stop("Not enough individuals in a centre.")
    
    dataset_centre <- dataset[indexes_centre,]
    dropout_matrix_centre <- dropout_matrix[indexes_centre,]
    e_matrix_centre <- e_matrix[indexes_centre,]
    
    # Compute the log-likelihood of the centre
    ll_centre <- ll_AdPaik_centre_1D(x, index, params, dataset_centre,
                                     dropout_matrix_centre, e_matrix_centre)
    ll_overall <- ll_overall + ll_centre
  }
  return (ll_overall)
}

#-------------------------------------------------------------------------------
#' @title
#' One-Dimensional Group log-Likelihood Function
#'
#' @description
#' This function simply implements the group log-likelihood function, following the definition.
#' It is internally used by @ll_AdPaik_1D and, therefore, it requires as first and second argument the parameter according to which the
#' global log-likelihood is one-dimensionally optimized and its position inside the vector of parameters.
#'
#' @param param_onedim One dimensional parameter, with respect to which the log-likelihood function
#' must be optmize.
#' @param index_param_onedim Index of the previous parameter inside the parameter vector.
#' @param params Parameter vector.
#' @param dataset Matrix of dataset regressors, with a number of rows equal to the number of individuals in a cluster.
#' @param dropout_matrix Binary matrix indicating in which interval of the time domain and individual failed. For an individual,
#' the sum of the row elements must be equal to 1 (if he/she failed) or 0 (if he/she does not failed).
#' It has dimension equal to (n_individuals, n_intervals)
#' @param e_matrix Matrix of dimension (n_individual, n_intervals), where each element contains the evaluation of the temporal
#' integral, performed through the function @param time_int_eval.
#'
#' @return Centre log-likelihood function.
#' 
#' @keywords internal


ll_AdPaik_centre_1D <- function(param_onedim, index_param_onedim, params, dataset,
                                dropout_matrix, e_matrix){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  R <- n_regressors <- dim(dataset)[2]                                       # Number of regressors
  L <- n_intervals <- dim(dropout_matrix)[2]                                 # Number of intervals
  n_params <- length(params)                                                 # Number of parameters
  
  # Impose value actual one_dim_parameter
  params[index_param_onedim] <- param_onedim
  
  # Extract parameters from the vector params
  phi <- matrix(params[1:L], nrow = L, ncol = 1)                                # Baseline log-hazard for L intervals
  betar <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)                      # Regression coefficients
  mu1 <- params[L+R+1]                                                          # Parameter for the gamma distribution of alpha
  nu <- params[L+R+2]                                                           # Parameter for the gamma distirbution of alpha
  gammak <- matrix(params[(L+3+R):(2*L+R+2)], nrow = L, ncol = 1)               # Parameter vector for the gamma distribution of epsilon
  
  # For idenfiability purposes
  mu2 <- 1 - mu1                                                                # Parameter for the gamma distribution of epsilon
  
  # Compute A_ijk
  A_ijk <- matrix(rep(0, n_individuals*L), nrow = n_individuals, ncol = L)
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% betar
    for (k in 1:L){
      A_ijk[j,k] <- as.numeric(e_matrix[j,k] * exp(data_betar + phi[k,1]))
    }
  }
  
  # Compute A_centre_dotdot
  #A_i <- colSums(matrix(rowSums(A_ijk), nrow = n_individuals, ncol = 1))
  A_i <- sum(A_ijk)
  
  # Compute A_ik
  A_ik <- matrix(colSums(A_ijk), nrow = 1, ncol = L)
  
  # Compute d_ik
  d_ik <- matrix(colSums(dropout_matrix), nrow = 1, ncol = L)
  
  # Compute the partial log-likelihood
  # First line of the formula
  loglik1 <- 0
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% betar
    for (k in 1:L){
      loglik1 <- loglik1 + as.numeric(dropout_matrix[j,k]*(data_betar + phi[k,1]))
    }
  }
  loglik1 <- loglik1 - as.numeric((mu1 / nu) * log(1 + nu * A_i))
  
  # Second line of the formula
  loglik2 <- 0
  for (k in 1:L){
    loglik2 <- loglik2 - as.numeric(mu2/gammak[k,1]) * log(1 + gammak[k,1] * A_ik[1,k])
  }
  
  # Third line of the formula
  loglik3 <- 0
  res_gamma1 <- gamma(mu1/nu)
  res1 <- (A_i + 1/nu)
  for (k in 1:L){
    loglik4 <- 0
    d_ikk <- d_ik[1,k]
    res_gamma2 <- gamma(mu2/gammak[k,1])
    res2 <- (d_ikk + mu2/gammak[k,1])
    res3 <- (A_ik[1,k] + 1/gammak[k,1])
    for (l in 0:d_ikk){
      coeff_binom <- choose(d_ikk, l)
      res_gamma3 <- gamma(res2 - l)
      res_gamma4 <- gamma(mu1/nu + l)
      res4 <- (res3)^(l - d_ikk)
      res5 <- (res1)^l
      res6 <- res_gamma4/res_gamma2
      if(res6 == 0)
        res6 <- 1e-10
      
      loglik4 <- loglik4 + res6 * (res_gamma3/res_gamma1) * (res4/res5) * coeff_binom
    }
    loglik3 <- loglik3 + log(loglik4)
  }
  
  # Return the sum of the three lines
  result <- loglik1 + loglik2 + loglik3
  return (result)
}

#-------------------------------------------------------------------------------
#' @title
#' Evaluation of Model log-Likelihood
#'
#' @description
#' Evaluation of the log-likelihood function at the provided parameter vector and data.
#'
#' @details
#' The function divides the individuals according to their group/cluster membership and then evaluates the group log-likelihood
#' through another implemented function, but using all and only the individuals belonging to that group.
#' Then the results are summed together to return the overall log-likelihood value.
#'
#' @param params Parameter vector
#' @param dataset Matrix of dimension equal to (number of individuals in the study, number of regressors), where only the regressors
#' indicated in the formula object are considered.
#' @param centre Vector of length equal to the number of individuals in the study, where each element corresponds to the individual
#' cluster membership.
#' @param time_axis Temporal domain
#' @param dropout_matrix Binary matrix indicating in which interval of the time domain and individual failed. For an individual,
#' the sum of the row elements must be equal to 1 (if he/she failed) or 0 (if he/she does not failed).
#' It has dimension equal to (n_individuals, n_intervals)
#' @param e_matrix Matrix of dimension (n_individual, n_intervals), where each element contains the evaluation of the temporal
#' integral, performed through the function @time_int_eval.
#'
#' @return Overall log-likelihood function value at the provided parameters and data
#'
#' @keywords internal

ll_AdPaik_eval <- function(params, dataset, centre, time_axis, dropout_matrix, e_matrix){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- length(params)
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Define a variable that contains the log-likelihood result
  # Then update it by summing the partial result of each centre (loglik_centre)
  ll_overall <- 0
  for (i in 1:n_centres){
    # Extract individuals in centre
    indexes_centre <- which(centre == centre_codes[i])
    dataset_centre <- dataset[indexes_centre,]
    e_matrix_centre <- e_matrix[indexes_centre,]
    dropout_matrix_centre <- dropout_matrix[indexes_centre,]
    
    # Compute the log-likelihood of the centre
    ll_centre <- ll_AdPaik_centre_eval(params, dataset_centre, dropout_matrix_centre, e_matrix_centre)
    ll_overall <- ll_overall + ll_centre
  }
  return (ll_overall)
}
#-------------------------------------------------------------------------------
#' @title
#' Evaluation of Model Group log-Likelihood
#'
#' @description
#' Evaluation of model group log-likelihood at the provided parameter vector and data.
#' This function is internally called by 'll_AdPaik_eval' to evaluate the log-likelihood function, considering
#' all and only the individuals belonging to a group.
#'
#' @param params Parameter vector.
#' @param dataset Matrix of dataset regressors, with a number of rows equal to the number of individuals in a cluster.
#' @param dropout_matrix Binary matrix indicating in which interval of the time domain and individual failed. For an individual,
#' the sum of the row elements must be equal to 1 (if he/she failed) or 0 (if he/she does not failed).
#' It has dimension equal to (n_individuals, n_intervals)
#' @param e_matrix Matrix of dimension (n_individual, n_intervals), where each element contains the evaluation of the temporal
#' integral, performed through the function @time_int_eval.
#'
#' @return Group log-likelihood evaluation
#'
#' @keywords internal

ll_AdPaik_centre_eval <- function(params, dataset, dropout_matrix, e_matrix){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  R <- n_regressors <- dim(dataset)[2]                                        # Number of regressors
  L <- n_intervals <- dim(dropout_matrix)[2]                                  # Number of intervals
  n_params <- length(params)                                                  # Number of parameters
  
  # Extract parameters from the vector params
  phi <- matrix(params[1:L], nrow = L, ncol = 1)                                # Baseline log-hazard for L intervals
  betar <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)                      # Regression coefficients
  mu1 <- params[L+R+1]                                                          # Parameter for the gamma distribution of alpha
  nu <- params[(L+R+2)]                                                         # Parameter for the gamma distirbution of alpha
  gammak <- matrix(params[(L+3+R):(2*L+R+2)], nrow = L, ncol = 1)               # Parameter vector for the gamma distribution of epsilon
  
  # For idenfiability purposes
  mu2 <- 1 - mu1                                                                # Parameter for the gamma distribution of epsilon
  
  # Compute
  # Compute A_ijk
  A_ijk <- matrix(rep(0, n_individuals*L), nrow = n_individuals, ncol = L)
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% betar
    for (k in 1:L){
      A_ijk[j,k] <- as.numeric(e_matrix[j,k] * exp(data_betar + phi[k,1]))
    }
  }
  
  # Compute A_centre_dotdot
  #A_i <- colSums(matrix(rowSums(A_ijk), nrow = n_individuals, ncol = 1))
  A_i <- sum(A_ijk)
  
  # Compute A_ik
  A_ik <- matrix(colSums(A_ijk), nrow = 1, ncol = L)
  
  # Compute d_ik
  d_ik <- matrix(colSums(dropout_matrix), nrow = 1, ncol = L)
  
  # Compute the partial log-likelihood
  # First line of the formula
  loglik1 <- 0
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% betar
    for (k in 1:L){
      loglik1 <- loglik1 + as.numeric(dropout_matrix[j,k]*(data_betar + phi[k,1]))
    }
  }
  loglik1 <- loglik1 - as.numeric((mu1 / nu) * log(1 + nu * A_i))
  
  # Second line of the formula
  loglik2 <- 0
  for (k in 1:L){
    loglik2 <- loglik2 - as.numeric(mu2/gammak[k,1]) * log(1 + gammak[k,1] * A_ik[1,k])
  }
  
  # Third line of the formula
  loglik3 <- 0
  res_gamma1 <- gamma(mu1/nu)
  res1 <- (A_i + 1/nu)
  for (k in 1:L){
    loglik4 <- 0
    d_ikk <- d_ik[1,k]
    res_gamma2 <- gamma(mu2/gammak[k,1])
    res2 <- (d_ikk + mu2/gammak[k,1])
    res3 <- (A_ik[1,k] + 1/gammak[k,1])
    for (l in 0:d_ikk){
      coeff_binom <- choose(d_ikk, l)
      res_gamma3 <- gamma(res2 - l)
      res_gamma4 <- gamma(mu1/nu + l)
      res4 <- (res3)^(l - d_ikk)
      res5 <- (res1)^l
      res6 <- res_gamma4/res_gamma2
      if(res6 == 0)
        res6 <- 1e-10
      
      loglik4 <- loglik4 + res6 * (res_gamma3/res_gamma1) * (res4/res5) * coeff_binom
    }
    loglik3 <- loglik3 + log(loglik4)
  }
  
  # Return the sum of the three lines
  result <- loglik1 + loglik2 + loglik3
  return (result)
}
