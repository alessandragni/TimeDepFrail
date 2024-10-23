StocTimeDepModel <- function(formula, data, time_axis,
                             categories_range_min, categories_range_max,
                             C_mult,
                             time_domain = 0, flag_time_domain = FALSE,
                             n_extrarun = 40, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                             print_previous_ll_values = c(TRUE, 3)){
  
  writeLines(sprintf("Stochastic Time-Dependent Centre-Specific Frailty Model:"))
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Check all input variables are provided
  if(missing(categories_range_max))
    stop("At least one input variable is missing, with no default.")
  
  # Check time_axis vector
  check.time_axis(time_axis)
  
  # Check elements of the dataset
  check.dataset(data)
  
  # Check positiveness of the multiplicative constant C
  check.C_mult(C_mult)
  
  # Check variables of formula are contained in dataset
  check.formula_terms(formula, data)
  
  # Extract elements of the formula object
  formula_variables <- all.vars(formula)
  
  # Extract position of cluster
  special <- c("cluster")
  terms_object <- terms(formula, special, data = data)
  cluster_index <- attr(terms_object, "specials")$cluster
  cluster_name <- formula_variables[cluster_index]
  
  # Extract response (time_to_event)
  response <- formula_variables[1]
  
  # Extract covariates
  covariates <- attr(terms_object, "term.labels")
  n_covariates <- length(covariates) - 1
  covariates <- covariates[1:n_covariates]
  
  # Create variables of the right structure: time_to_event, centre, dataset
  time_to_event <- as.numeric(data[,response])
  n_individuals <- length(time_to_event)
  
  centre <- data[,cluster_name]
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # In case of categorical variables
  dataset <- c()
  new_covariates <- c()
  for(j in 1:n_covariates){
    flag_char <- is.character(data[,covariates[j]])
    if(!flag_char){
      dataset <- cbind(dataset, as.numeric(data[,covariates[j]]))
      new_covariates <- c(new_covariates, covariates[j])
    }
    else{
      dummy_extracted <- extract_dummy_variables(data[,covariates[j]], covariates[j])
      dataset <- cbind(dataset, dummy_extracted$DummyMatrix)
      new_covariates <- c(new_covariates,dummy_extracted$DummyVariablesName)
    }
  }
  #dataset <- data.frame(dataset)
  
  # Extract information about n_regressors (R), n_intervals (L), n_individuals,
  # n_groups (N), n_params (P), n_run
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- n_intervals + n_regressors + 3
  n_run <- n_params + n_extrarun
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, 1, 1, 1)
  n_categories <- length(params_categories)
  
  # Check the correctness of provided range categories
  check.categories_params(n_categories, categories_range_min, categories_range_max)
  
  # Generate extended vector of parameters ranges
  params_range_min <- params_range_max <- c()
  for(c in 1: n_categories){
    n_params_in_c <- params_categories[c]
    params_range_min <- c(params_range_min, rep(categories_range_min[c], n_params_in_c))
    params_range_max <- c(params_range_max, rep(categories_range_max[c], n_params_in_c))
  }
  
  # Build the matrices e_{ijk} and d_{ijk}
  e_matrix <- matrix(rep(0, n_intervals * n_individuals), n_individuals, n_intervals)
  dropout_matrix <- matrix(rep(0, n_intervals * n_individuals), n_individuals, n_intervals)
  for(j in 1:n_individuals){
    for(k in 1:n_intervals){
      e_matrix[j,k] <- time_int_eval(time_to_event[j], k, time_axis)
      
      if ((time_to_event[j] < time_axis[k+1]) & (time_to_event[j] >= time_axis[k])){
        dropout_matrix[j,k] <- 1
      }
    }
  }
  
  # Initialize the vector of parameters
  params <- rep(0, n_params)
  for (p in 1:n_params){
    params[p] <- runif(1, params_range_min[p], params_range_max[p])
  }
  
  # Build the matrix containing the order according to which the log-likelihood function is
  # optimized at each run
  RunIndexes <- matrix(rep(0, n_run * n_params), n_run, n_params)
  for(i in 1:n_run){
    if(i <= n_params){    # Set the matrix to the ordered coefficient indexes
      actual_p <- (i-1)
      for(j in 1:n_params){
        actual_p <- actual_p + 1
        if(actual_p > n_params)
          actual_p <- 1
        
        RunIndexes[i,j] <- actual_p
      }
    }
    else{   # If exceed the number of parameters, set to casual values
      RunIndexes[i,] <- sample(seq(1, n_params), n_params, replace = FALSE)
    }
  }
  
  # Vector for the used and remaining indexes
  RemainingIndexes <- seq(1, n_params, 1)
  UsedIndexes <- c()
  
  # Matrix containing the optimized parameters and log-likelihood at each run
  global_optimal_params <- matrix(rep(0, n_run * n_params), n_run, n_params)
  global_optimal_loglikelihood <- rep(0, n_run)
  
  # Optimize the log-likelihood function
  writeLines(sprintf("Start log-likelihood optimization ... "))
  r <- 1                                                # Set the actual run
  actual_tol_ll <- 1                                    # Set the actual tolerance on the log-likelihood value
  ll_optimal <- -1e100                                  # Set initial value of the optimized log-likelihood to small value
  optimal_run <- 1                                      # Set initial value for optimal_run
  status <- TRUE                                        # Set TRUE to algorithm exit status
  
  while(r <= n_run & actual_tol_ll > tol_ll){
    writeLines(sprintf(paste("Run ", r)))
    
    # Select ordered indexes for current run
    RemainingIndexes <- RunIndexes[r,]
    UsedIndexes <- c()
    
    while(length(RemainingIndexes) != 0){
      # Select current index
      index_to_vary <- RemainingIndexes[1]
      PosIndex <- which(RemainingIndexes == index_to_vary)
      
      # Update remaining and used indexes
      RemainingIndexes <- RemainingIndexes[-PosIndex]
      UsedIndexes <- c(UsedIndexes,index_to_vary)
      
      result_optimize <- optimize(ll_StocTimeDep_1D,
                                  c(params_range_min[index_to_vary], params_range_max[index_to_vary]),
                                  maximum = TRUE, tol = tol_optimize,
                                  index_to_vary, params, dataset, centre,
                                  time_axis, dropout_matrix, e_matrix, time_to_event, C_mult)
      
      params[index_to_vary] <- result_optimize$maximum
    }
    
    global_optimal_params[r,] <- params
    global_optimal_loglikelihood_run <- ll_StocTimeDep_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix, time_to_event)
    global_optimal_loglikelihood[r] <- global_optimal_loglikelihood_run
    
    # Check meaningfulness of the global_optimal_loglikelihood
    if(is.nan(global_optimal_loglikelihood_run))
      stop("NaN value for the optimal log-likelihood value.")
    
    # Print previous values of the log-likelihood function
    if(print_previous_ll_values[1]){
      n_previous <- print_previous_ll_values[2]
      if(r < n_previous)
        writeLines(sprintf(paste("Global log-likelihood: ", global_optimal_loglikelihood[1:r])))
      else
        writeLines(sprintf(paste("Global log-likelihood: ", global_optimal_loglikelihood[(r - n_previous + 1):r])))
    }
    
    # Update conditions in while loop
    actual_tol_ll <- abs(ll_optimal - global_optimal_loglikelihood_run)
    if(ll_optimal < global_optimal_loglikelihood_run){
      ll_optimal <- global_optimal_loglikelihood_run
      optimal_run <- r
    }
    r <- r + 1
  }
  writeLines(sprintf(paste("... End optimization")))
  if(r == n_run)
    status = FALSE
  
  # Extract best solution with maximum log-likelihood
  optimal_params <- global_optimal_params[optimal_run,]
  optimal_loglikelihood <- global_optimal_loglikelihood[optimal_run]
  
  # Compute the standard error from the Hessian matrix
  writeLines(sprintf(paste("Compute parameters standard error")))
  params_se <- params_se.StocTimeDep(optimal_params, params_range_min, params_range_max,
                                     dataset, centre, time_axis, dropout_matrix, e_matrix, time_to_event, h_dd)
  
  # Compute parameters confidence interval
  writeLines(sprintf(paste("Compute parameters confidence interval")))
  params_CI <- params_CI(optimal_params, params_se)
  
  # Compute baseline hazard step-function
  writeLines(sprintf(paste("Compute baseline hazard step function")))
  bas_hazard <- bas_hazard(optimal_params, time_axis)
  
  # Compute frailty standard deviation
  writeLines(sprintf(paste("Compute frailty standard deviation")))
  check.time_domain(time_domain, flag_time_domain)
  if(flag_time_domain)
    frailty_dispersion <- frailty_sd.StocTimeDep(optimal_params, time_domain, n_regressors,
                                                 categories_range_min, categories_range_max)
  else
    frailty_dispersion <- frailty_sd.StocTimeDep(optimal_params, time_axis, n_regressors,
                                                 categories_range_min, categories_range_max)
  
  # Akaike Information Criterium
  AIC = 2 * n_param - 2 * optim_loglikelihood
  
  # Object to return
  return_list <- list("formula" = formula,
                      "Regressors" = new_covariates,
                      "NRegressors" = n_regressors,
                      "ClusterVariable" = cluster_name,
                      "NClusters" = n_centres,
                      "TimeDomain" = time_axis,
                      "NIntervals" = n_intervals,
                      "NParameters" = n_params,
                      "ParametersCategories" = params_categories,
                      "ParametersRange" = params_range,
                      "Loglikelihood" = optimal_loglikelihood,
                      "AIC" = AIC,
                      "Status" = status,
                      "NRun" = r-1,
                      "OptimalParameters" = optimal_params,
                      "StandardErrorParameters" = params_se,
                      "ParametersCI" = params_CI,
                      "BaselineHazard" = bas_hazard,
                      "FrailtyDispersion" = frailty_dispersion)
  class(return_list) <- "StocTimeDep"
  
  # Return list of results
  return (return_list)
}
#-------------------------------------------------------------------------------

# FUNCTIONS FOR IMPLEMETING THE LIKELIHOOD
ll_StocTimeDep_1D <- function(x, index, params, dataset, centre, time_axis,
                              dropout_matrix, e_matrix, time_to_event, C_mult){
  
  # Extract information from the input variables
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- length(params)
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Define a variable that contains the log-likelihood result.
  # It has an initial value corresponding to a constant quantity that does not depend on the centre.
  # Then update it by summing the partal result of each centre
  ll_overall <- (-n_centres * log(pi))
  for (i in 1:n_centres){
    # Extract individuals in a centre
    indexes_centre <- which(centre == centre_codes[i])
    
    # Check number of individuals in each centre
    if(length(indexes_centre) <= 1)
      stop("Not enough individuals in a centre.")
    
    dataset_centre <- dataset[indexes_centre,]
    dropout_matrix_centre <- dropout_matrix[indexes_centre,]
    e_matrix_centre <- e_matrix[indexes_centre,]
    time_to_event_centre <- time_to_event[indexes_centre]
    
    # Compute the log-likelihood of the centre
    ll_centre <- ll_StocTimeDep_centre_1D(x, index, params, dataset_centre, dropout_matrix_centre,
                                          e_matrix_centre, time_to_event_centre, C_mult)
    ll_overall <- ll_overall + ll_centre
  }
  return (ll_overall)
}
#-------------------------------------------------------------------------------
# Compute the log-likelihood referred to a centre
ll_StocTimeDep_centre_1D <- function(param_onedim, index_param_onedim, params, dataset,
                                     dropout_matrix, e_matrix, time_to_event, C_mult){
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  R <- n_regressors <- dim(dataset)[2]                                       # Number of regressors
  L <- n_intervals <- dim(dropout_matrix)[2]                                 # Number of intervals
  n_params <- length(params)                                                 # Number of parameters
  n_nodes <- dim(nodes_ghqm)[1]                                              # Number of nodes for the quadrature formula
  
  # Impose value actual one_dim_parameter
  params[index_param_onedim] <- param_onedim
  
  # Extract parameters from the vector params
  phi <- matrix(params[1:L], nrow = L, ncol = 1)                                # Baseline log-hazard for L intervals
  beta <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)                       # Regression coefficients
  lambda1 <- params[(L+R+1)]                                                    # Eigenvalues and angle alpha
  lambda2 <- params[(L+R+2)]
  angle_alpha <- params[(L+R+3)]
  
  # Get the variances and covariances
  sigma2c <- lambda1*(cos(angle_alpha))^2 + lambda2*(sin(angle_alpha))^2
  sigmacb <- (lambda1 - lambda2)*(sin(angle_alpha))*(cos(angle_alpha))
  sigma2b <- lambda1*(sin(angle_alpha))^2 + lambda2*(cos(angle_alpha))^2
  
  # Define further variables
  sigmac <- sqrt(sigma2c)
  sigmab <- sqrt(sigma2b)
  gamma <- sigmacb/(sigma2b)
  sigma2r <- sigma2c - sigma2b*(gamma^2)
  sigmar <- sqrt(sigma2r)
  
  # Compute
  # Compute d_ij
  d_ij <- matrix(rowSums(dropout_matrix), nrow = n_individuals, ncol = 1)
  
  # Compute d_i
  d_i <- colSums(d_ij)
  
  # Compute the partial log-likelihood
  # Compute first term of the formula
  loglik1 <- 0
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% beta
    for (k in 1:L){
      loglik1 <- loglik1 + as.numeric(dropout_matrix[j,k] * (data_betar + phi[k,1]))
    }
  }
  
  # Compute the third term of the formula
  loglik2 <- 0
  for(q in 1:n_nodes){
    arg_G <- nodes_ghqm[q,1]
    arg_exp <- sqrt(2) * sigmar * arg_G * d_i
    res1 <- weights_ghqm[q,1] * exp(arg_exp)
    G1 <- G_1D(param_onedim, index_param_onedim, params, dataset,
               time_to_event, time_axis, d_ij, d_i, arg_G)
    loglik2 <- loglik2 + as.numeric(res1* G1 * C_mult)
  }
  loglik2 <- log(loglik2)
  
  # Define the result
  result <- loglik1 + loglik2
  return (result)
}
#-------------------------------------------------------------------------------
# This function implements the G function of Wintrebert, pag 517 (correct version on the Thesis)
G_1D <- function(param_onedim, index_param_onedim, params, dataset,
                 time_to_event, time_axis, d_ij,d_i,z){
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Extract information
  n_individuals <- length(time_to_event)
  R <- n_regressors <- dim(dataset)[2]
  L <- n_intervals <- length(time_axis) - 1
  n_nodes <- dim(nodes_ghqm)[1]
  
  # Impose the single parameter
  params[index_param_onedim] <- param_onedim
  
  # Extract parameters from the vector
  phi <- matrix(params[1:L], nrow = L, ncol = 1)
  beta <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)
  lambda1 <- params[(L+R+1)]
  lambda2 <- params[(L+R+2)]
  angle_alpha <- params[(L+R+3)]
  
  # Get the variances and covariances
  sigma2c <- lambda1*(cos(angle_alpha))^2 + lambda2*(sin(angle_alpha))^2
  sigmacb <- (lambda1 - lambda2)*(sin(angle_alpha))*(cos(angle_alpha))
  sigma2b <- lambda1*(sin(angle_alpha))^2 + lambda2*(cos(angle_alpha))^2
  
  # Get the other variables
  sigmab <- sqrt(sigma2b)
  sigmac <- sqrt(sigma2c)
  gamma <- sigmacb/(sigma2b)
  sigma2r <- sigma2c - sigma2b*(gamma^2)
  sigmar <- sqrt(sigma2r)
  
  # Compute the function G(z)
  partial1 <- gamma*d_i
  partial2 <- sum(d_ij * time_to_event)
  partial <- 0
  for(u in 1:n_nodes){
    partial3 <- 0
    arg_exp1 <- sqrt(2)*sigmab*nodes_ghqm[u,1]
    for(jj in 1:n_individuals){
      data_betar <- as.numeric(dataset[jj,]) %*% beta
      time_to_event_j <- time_to_event[jj]
      for(kk in 1:L){
        res_f_ijk <- f_ijk(arg_exp1, phi, kk, time_to_event_j, time_axis)
        partial3 <- partial3 + exp(data_betar) * res_f_ijk
      }
    }
    arg_exp2 <- sqrt(2)*sigmar*z + arg_exp1*gamma
    arg_exp3 <- arg_exp1*(partial1+partial2) - partial3*(exp(arg_exp2))/arg_exp1
    partial <- partial + as.numeric(weights_ghqm[u,1]*exp(arg_exp3))
  }
  return (partial)
}
#-------------------------------------------------------------------------------
# This function implements the function f() of Wintrebert, pag 517
f_ijk <- function(b, phi, kkk, time_to_event_j, time_axis){
  if (time_to_event_j < time_axis[kkk])
    return (0)
  else if (time_to_event_j < time_axis[kkk+1] & time_to_event_j >= time_axis[kkk])
    return (exp(phi[kkk,1]) * (exp(b * time_to_event_j) - exp(b * time_axis[kkk])))
  else if(time_to_event_j >= time_axis[kkk+1])
    return (exp(phi[kkk,1]) * (exp(b * time_axis[kkk+1]) - exp(b * time_axis[kkk])))
}

#-------------------------------------------------------------------------------

# Compute the loglikelihood of the model
ll_StocTimeDep_eval <- function(params, dataset, centre, time_axis, dropout_matrix, e_matrix, time_to_event){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- length(params)
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Define a variable that contains the log-likelihood result.
  # It has an initial value corresponding to a constant quantity that does not depend on the centre.
  # Then update it by summing the partal result of each centre
  ll_overall <- (-n_centres * log(pi))
  for (i in 1:n_centres){
    # Extract individuals in centre
    indexes_centre <- which(centre == centre_codes[i])
    dataset_centre <- dataset[indexes_centre,]
    e_matrix_centre <- e_matrix[indexes_centre,]
    dropout_matrix_centre <- dropout_matrix[indexes_centre,]
    time_to_event_centre <- time_to_event[indexes_centre]
    
    # Compute the log-likelihood of the centre
    ll_centre <- ll_StocTimeDep_centre_eval(params, dataset_centre, dropout_matrix_centre,
                                            e_matrix_centre, time_to_event_centre)
    ll_overall <- ll_overall + ll_centre
  }
  return (ll_overall)
}
#-------------------------------------------------------------------------------

# Compute the log-likelihood referred to a centre
ll_StocTimeDep_centre_eval <- function(params, dataset, dropout_matrix, e_matrix,
                                       time_to_event){
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Extract information from input data
  n_individuals <- dim(dataset)[1]
  R <- n_regressors <- dim(dataset)[2]                                        # Number of regressors
  L <- n_intervals <- dim(dropout_matrix)[2]                                  # Number of intervals
  n_params <- length(params)                                                  # Number of parameters
  n_nodes <- dim(nodes_ghqm)[1]                                               # Number of nodes for the quadrature formula
  
  # Extract parameters from the vector params
  phi <- matrix(params[1:L], nrow = L, ncol = 1)                                # Baseline log-hazard for L intervals
  beta <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)                       # Regression coefficients
  lambda1 <- params[(L+R+1)]                                                    # Eigenvalues and angle alpha
  lambda2 <- params[(L+R+2)]
  angle_alpha <- params[(L+R+3)]
  
  # Get the variances and covariances
  sigma2c <- lambda1*(cos(angle_alpha))^2 + lambda2*(sin(angle_alpha))^2
  sigmacb <- (lambda1 - lambda2)*(sin(angle_alpha))*(cos(angle_alpha))
  sigma2b <- lambda1*(sin(angle_alpha))^2 + lambda2*(cos(angle_alpha))^2
  
  # Define further variables
  sigmac <- sqrt(sigma2c)
  sigmab <- sqrt(sigma2b)
  gamma <- sigmacb/(sigma2b)
  sigma2r <- sigma2c - sigma2b*(gamma^2)
  sigmar <- sqrt(sigma2r)
  
  # Compute
  # Compute d_ij
  d_ij <- matrix(rowSums(dropout_matrix), nrow = n_individuals, ncol = 1)
  
  # Compute d_i
  d_i <- colSums(d_ij)
  
  # Compute the partial log-likelihood
  # Compute first term of the formula
  loglik1 <- 0
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% beta
    for (k in 1:L){
      loglik1 <- loglik1 + as.numeric(dropout_matrix[j,k] * (data_betar + phi[k,1]))
    }
  }
  
  # Compute the third term of the formula
  loglik2 <- 0
  for(q in 1:n_nodes){
    arg_G <- nodes_ghqm[q,1]
    arg_exp <- sqrt(2) * sigmar * arg_G * d_i
    res1 <- weights_ghqm[q,1] * exp(arg_exp)
    G1_eval <- G_eval(params, dataset, time_to_event, time_axis,
                      d_ij, d_i, arg_G)
    loglik2 <- loglik2 + as.numeric(res1* G1_eval)
  }
  loglik2 <- log(loglik2)
  
  # Define the result
  result <- loglik1 + loglik2
  return (result)
}
#-------------------------------------------------------------------------------
# This function implements the G function of Wintrebert, pag 517 (correct version on the Thesis)
G_eval <- function(params, dataset, time_to_event, time_axis,
                   d_ij, d_i, z){
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Extract information from input variables
  n_individuals <- length(time_to_event)
  R <- n_regressor <- dim(dataset)[2]
  L <- n_intervals <- length(time_axis) - 1
  n_nodes <- dim(nodes_ghqm)[1]
  
  # Extract parameters
  phi <- matrix(params[1:L], nrow = L, ncol = 1)
  beta <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)
  lambda1 <- params[(L+R+1)]
  lambda2 <- params[(L+R+2)]
  angle_alpha <- params[(L+R+3)]
  
  # Get the variances and covariances
  sigma2c <- lambda1*(cos(angle_alpha))^2 + lambda2*(sin(angle_alpha))^2
  sigmacb <- (lambda1 - lambda2)*(sin(angle_alpha))*(cos(angle_alpha))
  sigma2b <- lambda1*(sin(angle_alpha))^2 + lambda2*(cos(angle_alpha))^2
  
  # Get the other variables
  sigmab <- sqrt(sigma2b)
  sigmac <- sqrt(sigma2c)
  gamma <- sigmacb/(sigma2b)
  sigma2r <- sigma2c - sigma2b*(gamma^2)
  sigmar <- sqrt(sigma2r)
  
  # Compute the function G(z)
  partial1 <- gamma*d_i
  partial2 <- sum(d_ij * time_to_event)
  partial <- 0
  for(u in 1:n_nodes){
    partial3 <- 0
    arg_exp1 <- sqrt(2)*sigmab*nodes_ghqm[u,1]
    for(jj in 1:n_individuals){
      data_betar <- as.numeric(dataset[jj,]) %*% beta
      time_to_event_j <- time_to_event[jj]
      for(kk in 1:L){
        res_f_ijk <- f_ijk(arg_exp1, phi, kk, time_to_event_j, time_axis)
        partial3 <- partial3 + exp(data_betar) * res_f_ijk
      }
    }
    arg_exp2 <- sqrt(2)*sigmar*z + arg_exp1*gamma
    arg_exp3 <- arg_exp1*(partial1+partial2) - partial3*(exp(arg_exp2))/arg_exp1
    partial <- partial + as.numeric(weights_ghqm[u,1]*exp(arg_exp3))
  }
  return (partial)
}
#-------------------------------------------------------------------------------
StocTimeDep_1D <- function(formula, data, time_axis,
                           index_param_to_vary, optimal_params = 0, flag_optimal_params = FALSE,
                           categories_range_min, categories_range_max,
                           C_mult,
                           n_iter = 5, tol_optimize = 1e-6,
                           flag_plot = TRUE, n_points = 150,
                           cex = 0.7, cex_max = 0.8, color_bg = "black", color_max_bg = "red",
                           pch = 21){
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Check all input variables are provided
  if(missing(C_mult))
    stop("At least one input variable is missing, with no default.")
  
  # Check positiveness of the multiplicative constant C
  check.C_mult(C_mult)
  
  # Check variables of formula are contained in dataset
  check.formula_terms(formula, data)
  
  # Check either the optimal parameters are provided or they must be simulated
  check.flag_optimal_params(optimal_params, flag_optimal_params)
  
  # Extract elements of the formula object
  formula_variables <- all.vars(formula)
  
  # Extract position of cluster
  special <- c("cluster")
  terms_object <- terms(formula, special, data = data)
  cluster_index <- attr(terms_object, "specials")$cluster
  cluster_name <- formula_variables[cluster_index]
  
  # Extract response (time_to_event)
  response <- formula_variables[1]
  
  # Extract covariates
  covariates <- attr(terms_object, "term.labels")
  n_covariates <- length(covariates) - 1
  covariates <- covariates[1:n_covariates]
  
  # Create variables of the right structure: time_to_event, centre, dataset
  time_to_event <- as.numeric(data[,response])
  n_individuals <- length(time_to_event)
  
  centre <- data[,cluster_name]
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # In case of categorical variables
  dataset <- c()
  new_covariates <- c()
  for(j in 1:n_covariates){
    flag_char <- is.character(data[,covariates[j]])
    if(!flag_char){
      dataset <- cbind(dataset, as.numeric(data[,covariates[j]]))
      new_covariates <- c(new_covariates, covariates[j])
    }
    else{
      dummy_extracted <- extract_dummy_variables(data[,covariates[j]], covariates[j])
      dataset <- cbind(dataset, dummy_extracted$DummyMatrix)
      new_covariates <- c(new_covariates,dummy_extracted$DummyVariablesName)
    }
  }
  #dataset <- data.frame(dataset)
  
  # Extract information about n_regressors (R), n_intervals (L), n_individuals,
  # n_groups (N), n_params (P)
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- n_intervals + n_regressors + 3
  
  # Check existence provided index
  check.index(index_param_to_vary, n_params)
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, 1, 1, 1)
  n_categories <- length(params_categories)
  
  # Check the correctness of provided range categories
  check.categories_params(n_categories, categories_range_min, categories_range_max)
  
  # Generate extended vector of parameters ranges
  params_range_min <- params_range_max <- c()
  for(c in 1: n_categories){
    n_params_in_c <- params_categories[c]
    params_range_min <- c(params_range_min, rep(categories_range_min[c], n_params_in_c))
    params_range_max <- c(params_range_max, rep(categories_range_max[c], n_params_in_c))
  }
  
  # Check that provided optimal parameters are contained in their min, max range
  if(flag_optimal_params)
    check.range_params(optimal_params, params_range_min, params_range_max)
  
  # Build the matrices e_{ijk} and d_{ijk}
  e_matrix <- matrix(rep(0, n_intervals * n_individuals), n_individuals, n_intervals)
  dropout_matrix <- matrix(rep(0, n_intervals * n_individuals), n_individuals, n_intervals)
  for(j in 1:n_individuals){
    for(k in 1:n_intervals){
      e_matrix[j,k] <- time_int_eval(time_to_event[j], k, time_axis)
      
      if ((time_to_event[j] < time_axis[k+1]) & (time_to_event[j] >= time_axis[k])){
        dropout_matrix[j,k] <- 1
      }
    }
  }
  
  # Store the estimated optimal parameter value and the optimal log-likelihood value
  param_optimal <- rep(0, n_iter)
  ll_optimized <- rep(0, n_iter)
  
  # Perform one-dimensional optimizatino
  params <- rep(0, n_params)
  for(iter in 1:n_iter){
    # Generate initial parameters according to the flag
    if(! flag_optimal_params){
      for(p in 1:n_params)
        params[p] <- runif(1, params_range_min[p], params_range_max[p])
    }
    else{
      params <- optimal_params
      params[index_param_to_vary] <- runif(1, params_range_min[index_param_to_vary], params_range_max[index_param_to_vary])
    }
    
    # Optimize the loglikelihood wrt indicated parameter
    result_optimize <- optimize(ll_StocTimeDep_1D,
                                c(params_range_min[index_param_to_vary], params_range_max[index_param_to_vary]),
                                maximum = TRUE, tol = tol_optimize,
                                index_param_to_vary, params, dataset, centre,
                                time_axis, dropout_matrix, e_matrix, time_to_event,
                                C_mult)
    
    param_optimal[iter] <- result_optimize$maximum
    ll_optimized[iter] <- result_optimize$objective
    
    if(flag_plot){
      param_1D <- param_optimal[iter]
      ll_1D <- ll_optimized[iter]
      plot_ll_1D.StocTimeDep(param_1D, index_param_to_vary, ll_1D, params,
                             params_range_min[index_param_to_vary], params_range_max[index_param_to_vary],
                             dataset, centre, time_axis, dropout_matrix, e_matrix, time_to_event,
                             n_points, cex, cex_max, color_bg, color_max_bg, pch)
    }
  }
  return_list <- list("EstimatedParameter" = param_optimal,
                      "OptimizedLoglikelihood" = ll_optimized)
  class(return_list) <- "StocTimeDep"
  
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
#' Function for plotting the trend of the log-likelihood function with respect to
#' a single parameter, specified in the arguments.
#' The log-likelihood function is plotted using a certain number of points specified by @n_points
#' and they are not connected by a straight line.
#'
#' @param_1D Optimal parameter value determined maximizing the log-likelihood function with respect to it.
#' @index_param_1D Index of the optimal parameter.
#' @ll_1D Log-likelihood value in correspondence of the optimal parameter, with other parameters fixed
#' to their optimal value or to a default value.
#' @params Value of the other parameters. A said, they could be fixed to a default value or to their optimal value.
#' @param_range_min Minimum value assumable by the parameter param_1D.
#' @param_range_max Maximum value assumable by the parameter param_1D.
#' @dataset Dataset with individual covariates. It can be either a matrix or a dataframe.
#' @centre Individual cluster membership. It is a vectorof the same dimension of the dataset.
#' @time_axis Temporal domain.
#' @dropout_matrix Binary matrix indicating in which interval of the time domain and individual failed. For an individual,
#' the sum of the row elements must be equal to 1 (if he/she failed) or 0 (if he/she does not failed).
#' It has dimension equal to (n_individuals, n_intervals)
#' @e_matrix Matrix of dimension (n_individual, n_intervals), where each element contains the evaluation of the temporal
#' integral, performed through the function @time_int_eval.
#' @n_points Number of points in which the log-likelihood function must be evaluated and then plotted.
#' To have a nice graphical representation, chose an intermediate value: not to small and not too high. Default value is 150.
#' @cex Dimension of the points (n_points). Deafult vaue is 0.7.
#' @cex_max Dimension of the optimal point. Default value is 0.8.
#' @color_bg Color of the points (n_points). Default is black.
#' @color_max_bg Color of the optimal point. Deafult is red.
#' @pch Shape of the points, with no distinction between optimal point and points.
plot_ll_1D.StocTimeDep <- function(param_1D, index_param_1D, ll_1D, params, param_range_min, param_range_max,
                                   dataset, centre, time_axis, dropout_matrix, e_matrix, time_to_event,
                                   n_points = 150,
                                   cex = 0.7, cex_max = 0.8, color_bg = "black", color_max_bg = "red",
                                   pch = 21){
  
  # Define the structure containing the generated points and the associated log-likelihood value
  param_values <- rep(0, n_points)
  ll_values <- rep(0, n_points)
  
  # Generate n_points for the indicated parameter inside its min, max range
  param_values <- runif(n_points, param_range_min, param_range_max)
  
  # For each point, evaluate the log-likelihood function
  for(i in 1:n_points){
    params[index_param_1D] <- param_values[i]
    ll_values[i] <- ll_StocTimeDep_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix, time_to_event)
  }
  
  # Plot the log-likelihood trend with respect to the indicated parameter
  string_title <- paste("Log-likelihood trend wrt parameter ", index_param_1D)
  
  # dev.new()
  plot(param_values, ll_values, pch=pch, col=color_bg, cex = cex,
       xlim = c(param_range_min, param_range_max), ylim=c(min(ll_values), max(ll_values)),
       main = string_title, xlab = "Values", ylab = "Log-likelihood")
  points(param_1D, ll_1D, bg = color_max_bg, pch = pch, cex = cex_max)
}


#-------------------------------------------------------------------------------
#' @title
#' Summary of the 'Adapted Paik et al.'s Model'
#'
#' @description
#' Summary function for summarizing the most important information related to the dataset (number of individuals,
#' number of regressors, number of intervals, number of clusters), the model call (number of parameters) and the model
#' output (optimal log-likelihood value and AIC).
#'
#' @details
#' Among the estimated parameters, only the regressors are reported together with their standard error and confidence interval.
#'
#' @param result 'S3' class object returned by the main model call, i.e. output of the 'Adapted Paik et al.'s Model'.
#'
#' @return Model summary printed on output.
#'
#' @export
#'
#' @examples
#' # Define the variables needed for the model execution
#' formula <-
#' time_axis <-
#' categories_range_min <- c()
#' categories_range_max <- c()
#'
#' # Call the main model function
#' result <- StocTimeDep(formula, data, time_axis, categories_range_min, categories_range_max)
#'
#' # Call the summary
#' summary.StocTimeDep(result)

summary.StocTimeDep <- function(result){
  check.result(result)
  
  # Extract information from the model output
  params_categories <- result$ParametersCategories
  n_categories <- length(params_categories)
  L <- n_intervals <- params_categories[1]
  R <- n_regressors <- params_categories[2]
  
  # Create new vector where each optimal parameter is followed by its standard error
  n_params <- result$NParameters
  optimal_parameters <- rep(0, n_params)
  for(p in 1:n_params){
    optimal_parameters[p] <- paste(round(result$OptimalParameters[p],4), round(result$StandardErrorParameters[p],4), sep = " (")
    optimal_parameters[p] <- paste(optimal_parameters[p], "", sep=")")
  }
  
  # Initialize vector for estimated regressors
  betar <- optimal_parameters[(L+1):(L+R)]
  
  # Create other variables for the ouput
  convergence <- ""
  if(result$Status == TRUE){
    convergence <- paste("TRUE (Convergence in ", result$NRun)
    convergence <- paste(convergence, " runs).")
  }else
    convergence <- "FALSE (No Convergence)"
  
  string_parameters <- paste("Overall number of parameters ", result$NParameters)
  string_parameters <- paste(string_parameters, "divided as (phi, betar, sigmac, sigmab, sigmacb) = (", sep=",\n")
  for(p in 1:n_categories){
    if(p == n_categories)
      string_parameters <- paste(string_parameters, params_categories[p],")")
    else
      string_parameters <- paste(string_parameters, params_categories[p],",")
  }
  
  # Extract entire formula call
  formula_string <- paste(result$formula[2], result$formula[1], result$formula[3])
  
  # Print output
  paste0 <- paste("Call: ", formula_string)
  paste9 <- paste("with cluster variable '",result$ClusterVariable,"' (", result$NClusters,"clusters).")
  paste1 <- paste("Log-likelihood:           ", round(result$Loglikelihood,4))
  paste2 <- paste("AIC:                       ", round(result$AIC,4))
  paste3 <- paste("Status of the algorithm:   ", convergence)
  paste4 <- "-------------------------------------------------------------------------------"
  paste5 <- string_parameters
  paste6 <- paste("with: number of intervals =", n_intervals)
  paste7 <- paste("      number of regressors =", n_regressors, ".")
  paste8 <- paste("Estimated regressors (standard error):")
  
  output <- paste("Output of the 'Stochastic Time-Dependent Centre-Specific Frailty Model'", paste4, sep="\n")
  output <- paste(output, paste0, sep="\n")
  output <- paste(output, paste9, sep="\n")
  output <- paste(output, paste4, sep="\n")
  #--------------
  output <- paste(output, paste1, sep="\n")
  output <- paste(output, paste2, sep="\n")
  output <- paste(output, paste3, sep="\n")
  output <- paste(output, paste4, sep="\n")
  #--------------
  output <- paste(output, paste5, sep="\n")
  output <- paste(output, paste6, sep=",\n")
  output <- paste(output, paste7, sep="\n")
  output <- paste(output, paste4, sep="\n")
  #-------------
  output <- paste(output, paste8, sep="\n")
  for(r in 1:R){
    string_regressor <- paste(result$Regressors[r],":",betar[r])
    output <- paste(output, string_regressor, sep="\n")
  }
  output <- paste(output, paste4, sep="\n")
  cat(output)
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
check.result.StocTimeDep <- function(result){
  # Save the names of the list elements
  names_list.StocTimeDep <- c("formula", "Regressors", "NRegressors", "ClusterVariable", "NClusters",
                              "TimeDomain", "NIntervals",
                              "NParameters", "ParametersCategories",
                              "ParametersRange",
                              "Loglikelihood", "AIC", "Status", "NRun",
                              "OptimalParameters", "StandardErrorParameters",
                              "ParametersCI", "BaselineHazard",
                              "FrailtyDispersion")
  
  # Other than a class, it is a list
  if(! is.list(result))
    stop("Wrong structure for input 'result' argument.")
  
  names_list <- names_list.StocTimeDep
  for(i in 1:length(names_list)){
    if(names(result)[i] != names_list[i]){
      msg <- paste(names_list[i], "does not appear in the input 'result' argument. ")
      stop(msg)
    }
  }
  
  # Compute the number of parameters
  n_params <- result$NIntervals + result$NRegressors + 3
  
  # For each element of the list, control its structure
  if(class(result$formula) != "formula")
    stop("'formula' is not a formula object.")
  if(!is.vector(result$Regressors))
    stop("'Regressors' is not a vector.")
  if(!is.numeric(result$NRegressors))
    stop("'NRegressors' is not a number.")
  if(!is.character(result$ClusterVariable))
    stop("'ClusterVariable' is not a string.")
  if(!is.numeric(result$NClusters))
    stop("'NCluster' is not a number.")
  
  if(! is.vector(result$TimeDomain))
    stop("'TimeDomain' is not a vector.")
  if(!is.numeric(result$NIntervals))
    stop("'NIntervals' is not a number.")
  if(length(result$TimeDomain) - 1 != result$NIntervals)
    stop("Different values for number of intervals in 'TimeDomain' and 'NIntervals'")
  
  if(! is.numeric(result$NParameters))
    stop("'NParameters' is not a number.")
  if(! is.vector(result$ParametersCategories))
    stop("'ParametersCategories' is not a vector.")
  if(length(result$ParametersCategories) != 5)
    stop("Wrong length of 'ParametersCategories' vector.")
  
  check.params_range(result$ParametersRange, n_params)
  
  if(! is.numeric(result$Loglikelihood))
    stop("'Loglikelihood' is not a value.")
  if(! is.numeric(result$AIC))
    stop("'AIC' is not a value.")
  if(! is.logical(result$Status))
    stop("'Status' is not a binary variable.")
  if(! is.numeric(result$NRun))
    stop("NRun' is not a number.")
  
  if(! is.vector(result$OptimalParameters))
    stop("'OptimalParameters' is not a vector.")
  if(! is.vector(result$StandardErrorParameters))
    stop("'StandardErrorParameters' is not a vector.")
  if(length(result$OptimalParameters) != n_params)
    stop("Wrong length of 'OptimalParameters' vector.")
  if(length(result$StandardErrorParameters) != n_params)
    stop("Wrong length of 'StandardErrorParameters' vector.")
  
  check.structure_paramsCI(result$ParametersCI)
  
  if(! is.vector(result$BaselineHazard))
    stop("'BaselineHazard' is not a vector.")
  if(length(result$BaselineHazard) != result$NIntervals)
    stop("Wrong length of 'BaselineHazrad' vector.")
  
  check.frailty_dispersion(result$FrailtyDispersion, result$NIntervals)
}


#-------------------------------------------------------------------------------
params_se.StocTimeDep <- function(optimal_params, params_range_min, params_range_max,
                                  dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd){
  
  # Assign nodes and weights
  nodes_ghqm <- nodesG_ghqm
  weights_ghqm <- weightsG_ghqm
  
  # Extract information from input variables
  n_params <- length(optimal_params)
  n_intervals <- length(time_axis) - 1
  n_regressors <- dim(dataset)[2]
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, 1, 1, 1)
  n_categories <- length(params_categories)
  
  # Generate extended vector of parameters ranges
  params_range_min <- params_range_max <- c()
  for(c in 1: n_categories){
    n_params_in_c <- params_categories[c]
    params_range_min <- c(params_range_min, rep(categories_range_min[c], n_params_in_c))
    params_range_max <- c(params_range_max, rep(categories_range_max[c], n_params_in_c))
  }
  
  # Initialize parameters standard error vector
  se <- rep(0, n_params)
  
  # Initialize information and hessian element
  information_element <- hessian_element <- 0
  
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
    
    ll_eval <- ll_StocTimeDep_eval(optimal_params, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_eval_plus <- ll_StocTimeDep_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_eval_minus <- ll_StocTimeDep_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    
    # Approximate the second derivative of the log-likelihood function
    hessian_element <- (ll_eval_plus + ll_eval_minus - 2*ll_eval)/(h_dd * h_dd)
    
    if((hessian_element == Inf) || (hessian_element == -Inf)){
      #information_element <- hessian_element
      se[p] <- 1e-4
    }
    else if(is.nan(hessian_element)){
      se[p] <- NaN
    }
    else{
      # Compute the information element from the hessian
      information_element <- - hessian_element
      
      # Compute standard error of the parameter
      se[p] <- 1/sqrt(information_element)
    }
  }
  return (se)
}
