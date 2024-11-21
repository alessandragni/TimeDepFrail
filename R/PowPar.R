PowParModel <- function(formula, data, time_axis,
                        categories_range_min, categories_range_max,
                        C_mult,
                        n_extrarun = 60, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                        print_previous_ll_values = c(TRUE, 3),
                        verbose = FALSE){
  if (verbose) message("Centre-Specific Frailty Model with Power Parameter:")
  
  # Assign nodes and weights
  nodes_ghqm <- nodes9_ghqm
  weights_ghqm <- weights9_ghqm
  
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
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- 2 * n_intervals + n_regressors
  #n_run <- n_params + n_extrarun
  n_run <- 1  # just to try it
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, n_intervals - 1, 1)
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
  params_range <- list("ParametersRangeMin" = params_range_min,
                       "ParametersRangeMax" = params_range_max)
  class(params_range) <- "ParametersRange"
  
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
  if (verbose) message("Start log-likelihood optimization ... ")
  r <- 1                                                # Set the actual run
  actual_tol_ll <- 1                                    # Set the actual tolerance on the log-likelihood value
  ll_optimal <- -1e100                                  # Set initial value of the optimized log-likelihood to small value
  optimal_run <- 1                                      # Set initial value for optimal_run
  status <- TRUE                                        # Set TRUE to algorithm exit status
  
  while(r <= n_run & actual_tol_ll > tol_ll){
    if (verbose) message(paste("Run ", r))
    
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
      
      # Optimize the log-likelihood wrt indicated index index_to_vary
      result_optimize <- optimize(ll_PowPar_1D,
                                  c(params_range_min[index_to_vary], params_range_max[index_to_vary]),
                                  maximum = TRUE, tol = tol_optimize,
                                  index_to_vary, params, dataset, centre,
                                  time_axis, dropout_matrix, e_matrix,
                                  C_mult)
      
      params[index_to_vary] <- result_optimize$maximum
    }
    
    global_optimal_params[r,] <- params
    global_optimal_loglikelihood_run <- ll_PowPar_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix)
    global_optimal_loglikelihood[r] <- global_optimal_loglikelihood_run
    
    # Check meaningfulness of the global_optimal_loglikelihood
    if(is.nan(global_optimal_loglikelihood_run))
      stop("NaN value for the optimal log-likelihood value.")
    
    # Print previous values of the log-likelihood function
    if(print_previous_ll_values[1]){
      n_previous <- print_previous_ll_values[2]
      if(r < n_previous)
        if (verbose) message(paste("Global log-likelihood: ", global_optimal_loglikelihood[1:r]))
      else
        if (verbose) message(paste("Global log-likelihood: ", global_optimal_loglikelihood[(r - n_previous + 1):r]))
    }
    
    # Update conditions in while loop
    actual_tol_ll <- abs(ll_optimal - global_optimal_loglikelihood_run)
    if(ll_optimal < global_optimal_loglikelihood_run){
      ll_optimal <- global_optimal_loglikelihood_run
      optimal_run <- r
    }
    r <- r + 1
  }
  if (verbose) message(paste("... End optimization"))
  if(r == n_run)
    status = FALSE
  
  # Extract best solution with maximum log-likelihood
  optimal_params <- global_optimal_params[optimal_run,]
  optimal_loglikelihood <- global_optimal_loglikelihood[optimal_run]
  
  # Compute the standard error from the Hessian matrix
  if (verbose) message(paste("Compute parameters standard error"))
  params_se <- params_se.PowPar(optimal_params, params_range_min, params_range_max,
                                dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd)
  
  # Compute parameters confidence interval
  if (verbose) message(paste("Compute parameters confidence interval"))
  params_CI <- params_CI(optimal_params, params_se)
  
  # Compute baseline hazard step-function
  if (verbose) message(paste("Compute baseline hazard step function"))
  bas_hazard <- bas_hazard(optimal_params, time_axis)
  
  # Compute frailty standard deviation
  if (verbose) message(paste("Compute frailty standard deviation"))
  frailty_dispersion <- frailty_sd.PowPar(optimal_params, time_axis, n_regressors,
                                          categories_range_min, categories_range_max)
  
  # Compute estimated frailty mean
  if (verbose) message(paste("Compute estimated frailty mean"))
  frailty_mean <- frailty_mean.PowPar(optimal_params, time_axis, n_regressors)
  
  # Akaike Information Criterium
  AIC = 2 * n_params - 2 * optimal_loglikelihood
  
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
                      "FrailtyDispersion" = frailty_dispersion,
                      "EstimatedFrailtyMean" = frailty_mean)
  class(return_list) <- "PowPar"
  
  # Return list of results
  return (return_list)
}

#------------------------------------------------------------------------------
ll_PowPar_1D <- function(x, index, params, dataset, centre, time_axis,
                         dropout_matrix, e_matrix, C_mult){
  
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
    ll_centre <- ll_PowPar_centre_1D(x, index, params, dataset_centre, dropout_matrix_centre,
                                     e_matrix_centre, C_mult)
    ll_overall <- ll_overall + ll_centre
  }
  
  # Return the resulted loglikelihood, but removing the constant term (second term of first line of the formula)
  return (ll_overall -(n_centres/2)*log(pi))
}

#-------------------------------------------------------------------------------
# Compute the log-likelihood referred to a centre
ll_PowPar_centre_1D <- function(param_onedim, index_param_onedim, params, dataset,
                                dropout_matrix, e_matrix, C_mult){
  
  # Assign nodes and weights
  nodes_ghqm <- nodes9_ghqm
  weights_ghqm <- weights9_ghqm
  
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
  betar <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)                      # regression coefficients
  gammak <- matrix(rep(0,L), nrow = L, ncol = 1)                                # frailty time-dependence parameter.
  gammak[2:L,1] <- params[(L+1+R):(2*L+R-1)]
  sigma <- params[2*L+R]                                                        # standard deviation of normal Yi
  
  # For idenfiability purposes
  gammak[1,1] <- 1
  
  # Compute log-likelihood:
  # Compute the first term of the first line of the formula
  loglik1 <- 0
  partial1 <- 0
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% betar
    for (k in 1:L){
      loglik1 <- loglik1 + as.numeric(dropout_matrix[j,k] * (data_betar + phi[k,1]))
      partial1 <- partial1 + as.numeric(dropout_matrix[j,k]*gammak[k,1])
    }
  }
  
  # Compute the second line
  loglik2 <- 0
  for(q in 1:n_nodes){
    partial2 <- 0
    for(j in 1:n_individuals){
      data_betar <- as.numeric(dataset[j,]) %*% betar
      for(k in 1:L){
        partial2 <- (partial2 + as.numeric(e_matrix[j,k]*exp(sqrt(2)*sigma*gammak[k,1]*nodes_ghqm[q,1] + phi[k,1] + data_betar)))
      }
    }
    loglik2 <- loglik2 + as.numeric(weights_ghqm[q,1]*exp(sqrt(2)*sigma*nodes_ghqm[q,1]*partial1-partial2))*(C_mult)
  }
  loglik2 <- log(loglik2)
  result <- loglik1 + loglik2
  return (result)
}

#-------------------------------------------------------------------------------

ll_PowPar_eval <- function(params, dataset, centre, time_axis, dropout_matrix, e_matrix){
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  n_regressors <- dim(dataset)[2]
  n_intervals <- length(time_axis) - 1
  n_params <- length(params)
  
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # Define a variable that contains the log-likelihood result
  # Then update it by summing the partial result of each centre
  ll_overall <- 0
  for (i in 1:n_centres){
    # Extract individuals in centre
    indexes_centre <- which(centre == centre_codes[i])
    dataset_centre <- dataset[indexes_centre,]
    e_matrix_centre <- e_matrix[indexes_centre,]
    dropout_matrix_centre <- dropout_matrix[indexes_centre,]
    
    # Compute the log-likelihood of the centre
    ll_centre <- ll_PowPar_centre_eval(params, dataset_centre, dropout_matrix_centre, e_matrix_centre)
    ll_overall <- ll_overall + ll_centre
  }
  
  # Return the resulted loglikelihood, but removing the constant term (second term of first line of the formula)
  return (ll_overall -(n_centres/2)*log(pi))
}
#-------------------------------------------------------------------------------

# Compute the log-likelihood referred to a centre
ll_PowPar_centre_eval <- function(params, dataset, dropout_matrix, e_matrix){
  
  # Assign nodes and weights
  nodes_ghqm <- nodes9_ghqm
  weights_ghqm <- weights9_ghqm
  
  # Extract information from input variables
  n_individuals <- dim(dataset)[1]
  R <- n_regressors <- dim(dataset)[2]                                        # Number of regressors
  L <- n_intervals <- dim(dropout_matrix)[2]                                  # Number of intervals
  n_params <- length(params)                                                  # Number of parameters
  n_nodes <- dim(nodes_ghqm)[1]                                               # Number of nodes for the quadrature formula
  
  # Extract parameters from the vector params
  phi <- matrix(params[1:L], nrow = L, ncol = 1)                                # Baseline log-hazard for L intervals
  betar <- matrix(params[(L+1):(L+R)], nrow = R, ncol = 1)                      # Regression coefficients
  gammak <- matrix(rep(0,L), nrow = L, ncol = 1)                                # Frailty time-dependence parameter.
  gammak[2:L,1] <- params[(L+1+R):(2*L+R-1)]
  sigma <- params[2*L+R]                                                        # standard deviation of normal Yi
  
  # For idenfiability purposes
  gammak[1,1] <- 1
  
  # Compute log-likelihood:
  # Compute the first term of the first line of the formula
  loglik1 <- 0
  partial1 <- 0
  for (j in 1:n_individuals){
    data_betar <- as.numeric(dataset[j,]) %*% betar
    for (k in 1:L){
      loglik1 <- loglik1 + as.numeric(dropout_matrix[j,k] * (data_betar + phi[k,1]))
      partial1 <- partial1 + as.numeric(dropout_matrix[j,k]*gammak[k,1])
    }
  }
  
  # Compute the second line
  loglik2 <- 0
  for(q in 1:n_nodes){
    partial2 <- 0
    for(j in 1:n_individuals){
      data_betar <- as.numeric(dataset[j,]) %*% betar
      for(k in 1:L){
        partial2 <- (partial2 + as.numeric(e_matrix[j,k]*exp(sqrt(2)*sigma*gammak[k,1]*nodes_ghqm[q,1] + phi[k,1] + data_betar)))
      }
    }
    loglik2 <- loglik2 + as.numeric(weights_ghqm[q,1]*exp(sqrt(2)*sigma*nodes_ghqm[q,1]*partial1-partial2))
  }
  loglik2 <- log(loglik2)
  result <- loglik1 + loglik2
  
  return (result)
}
#-------------------------------------------------------------------------------

# CSFM WITH POWER PARAMETER
PowPar_1D <- function(formula, data, time_axis,
                      index_param_to_vary, optimal_params = 0, flag_optimal_params = FALSE,
                      categories_range_min, categories_range_max,
                      C_mult,
                      n_iter = 5, tol_optimize = 1e-6,
                      flag_plot = TRUE, n_points = 150,
                      cex = 0.7, cex_max = 0.8, color_bg = "black", color_max_bg = "red",
                      pch = 21){
  
  # Assign nodes and weights
  nodes_ghqm <- nodes9_ghqm
  weights_ghqm <- weights9_ghqm
  
  # Check all input variables are provided
  if(missing(C_mult))
    stop("At least one input variable is missing, with no default.")
  
  # Check variables of formula are contained in dataset
  check.formula_terms(formula, data)
  
  # Check positiveness of the multiplicative constant C
  check.C_mult(C_mult)
  
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
  n_params <- 2 * n_intervals + n_regressors
  
  # Check existence provided index
  check.index(index_param_to_vary, n_params)
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, n_intervals - 1, 1)
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
  param_to_optimize <- rep(0, n_iter)
  ll_to_be_optimized <- rep(0, n_iter)
  
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
    result_optimize <- optimize(ll_PowPar_1D,
                                c(params_range_min[index_param_to_vary], params_range_max[index_param_to_vary]),
                                maximum = TRUE, tol = tol_optimize,
                                index_param_to_vary, params, dataset, centre,
                                time_axis, dropout_matrix, e_matrix,
                                C_mult)
    
    param_to_optimize[iter] <- result_optimize$maximum
    ll_to_be_optimized[iter] <- result_optimize$objective
    
    if(flag_plot){
      param_1D <- param_to_optimize[iter]
      ll_1D <- ll_to_be_optimized[iter]
      plot_ll_1D.PowPar(param_1D, index_param_to_vary, ll_1D, params,
                        params_range_min[index_param_to_vary], params_range_max[index_param_to_vary],
                        dataset, centre, time_axis, dropout_matrix, e_matrix,
                        n_points, cex, cex_max, color_bg, color_max_bg, pch)
    }
  }
  return_list <- list("EstimatedParameter" = param_optimal,
                      "OptimizedLoglikelihood" = ll_optimized)
  class(return_list) <- "PowPar_1D"
  
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
plot_ll_1D.PowPar <- function(param_1D, index_param_1D, ll_1D, params, param_range_min, param_range_max,
                              dataset, centre, time_axis, dropout_matrix, e_matrix,
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
    ll_values[i] <- ll_PowPar_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix)
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
summary.PowPar <- function(result){
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
  string_parameters <- paste(string_parameters, "divided as (phi, betar, gammak, sigma) = (", sep=",\n")
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
  
  output <- paste("Output of the 'Centre-Specific Frailty Model with Power Parameter'", paste4, sep="\n")
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
check.result.PowPar <- function(result){
  # Save the names of the list elements
  names_list.PowPar <- c("formula", "Regressors", "NRegressors", "ClusterVariable", "NClusters",
                         "TimeDomain", "NIntervals",
                         "NParameters", "ParametersCategories",
                         "ParametersRange",
                         "Loglikelihood", "AIC", "Status", "NRun",
                         "OptimalParameters", "StandardErrorParameters",
                         "ParametersCI", "BaselineHazard",
                         "FrailtyDispersion", "EstimatedFrailtyMean")
  
  # Other than a class, it is a list
  if(! is.list(result))
    stop("Wrong structure for input 'result' argument.")
  
  names_list <- names_list.PowPar
  for(i in 1:length(names_list)){
    if(names(result)[i] != names_list[i]){
      msg <- paste(names_list[i], "does not appear in the input 'result' argument. ")
      stop(msg)
    }
  }
  
  # Compute the number of parameters
  n_params <- result$NIntervals * 2 + result$NRegressors
  
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
  if(length(result$ParametersCategories) != 4)
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
  
  if (! is.vector(result$EstimatedFrailtyMean))
    stop(" 'EstimatedFrailtyMean' is not a vector.")
  if(length(result$EstimatedFrailtyMean) != result$NIntervals)
    stop("Wrong length of 'EstimatedFrailtyMean' vector.")
}


#-------------------------------------------------------------------------------

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


#-------------------------------------------------------------------------------
# Functions for computing the second derivative
params_se.PowPar <- function(optimal_params, params_range_min, params_range_max,
                             dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd){
  
  # Assign nodes and weights
  nodes_ghqm <- nodes9_ghqm
  weights_ghqm <- weights9_ghqm
  
  # Extract information from input variables
  n_params <- length(optimal_params)
  n_intervals <- length(time_axis) - 1
  n_regressors <- dim(dataset)[2]
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, n_intervals - 1, 1)
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
    
    ll_eval <- ll_PowPar_eval(optimal_params, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_eval_plus <- ll_PowPar_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_eval_minus <- ll_PowPar_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    
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

