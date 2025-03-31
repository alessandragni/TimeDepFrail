#' @title
#' Adapted Paik et Al.'s Model: Time-Dependent Shared Frailty Cox Model
#'
#' @description
#' Function for applying the 'Adapted Paik et al.'s Model', a Cox Model with time-dependent frailty,
#' shared by individuals belonging to the same group/cluster.
#'
#' To generate time-dependence, the temporal domain is divided into a certain number
#' of intervals. The model log-likelihood function depends on a certain number of parameters and 
#' is maximized with respect to all of them,
#' using a reinterpretation of the 'Powell's method in multidimension', 
#' that is a multi-dimensional optimization method based on
#' repeated one-dimensional optimization of the log-likelihood function (with respect to one parameter at the time).
#' In this context, the one-dimensional optimization is performed through the 'optimize' R function.
#' For more information about the unknown model parameters, their type and numerosity refer to Details.
#'
#' Several quantities are estimated at the end of the optimization phase, such as
#' optimal parameters, baseline hazard, frailty dispersion (standard deviation and variance),
#' posterior frailty estimates, with their variance and confidence interval, 
#' conditional survival function, Akaike Information Criterion (AIC), ...
#'
#' @details
#' Two observation needs to made about the time-domain:
#' - The time domain may coincide with the follow-up or it may be contained in it.
#' Indeed, the left boundary can be greater than the beginning of the follow-up and, for instance, it can coincide
#' with the time-instants in which the events begin to happen; conversely, the right boundary of the two must be the same.
#' - The partition of the time domain into intervals can be made according to two selected criteria:
#' (1) using an already existent partition of the follow-up
#' (2) using the shape of the baseline hazard function for a time independent model as reference: divide the time-domain according to regions in
#' which it has a peak or a plateau.
#'
#' @details
#' The parameters with respect to which the log-likelihood function must be optimized are:
#' - baseline log-hazard (number of parameters = number of intervals of the time-domain)
#' - data regressors
#' - \eqn{\mu_1}, \eqn{\nu}: parameters of the gamma distribution of \eqn{\alpha_j} (time-independent/constant) (2 parameters)
#' - \eqn{\gamma_k}: parameters of the gamma distribution of \eqn{\epsilon_{jk}} (time-dependent) (number of parameters = number of intervals)
#' Another model parameter is \eqn{\mu_2} and it is get imposing the constraint that \eqn{\mu_1 + \mu_2 = 1}.
#' As it can be notice, some parameters can be grouped into the same category (regressors, baseline log-hazard and so on)
#' and we can easily constraint them assigning each category both a minimum and maximum range.
#' The vector is structured as follows: (baseline log-hazard, regressors, \eqn{\mu_1}, \eqn{\nu}, \eqn{\gamma_k}) with dimension
#' (n_intervals, n_regressors, 1, 1, n_intervals).
#'
#'
#' @param formula Formula object having on the left hand side the time-to-event variable, that is the time-instant in which
#' the individual failed. On the right hand side, it has the regressors and the cluster variable.
#' @param data Dataset in which all variables of the formula object must be found and contained.
#' This dataset can also contain other variables, but they will not be considered.
#' It can be either a dataframe or a matrix, but in both cases the name of each column must correspond to the name of
#' the formula variables. It is not necessary to attach it (in case of a data.frame)
#' @param time_axis Temporal domain
#' @param categories_range_min Vector containing the minimum value assumable by each parameter category.
#' @param categories_range_max Vector containing the maximum value assumable by each parameter category.
#' @param n_extrarun Total number of runs (iterations) are obtained summing to the number of parameters and n_extrarun.
#' @param tol_ll Tolerance on the log-likelihood value.
#' @param tol_optimize Internal tolerance for the one-dimensional optimization through 'optimize' R function.
#' @param h_dd Discretization step used for the numerical approximation of the second derivative of the log-likelihood function.
#' @param verbose Logical. If `TRUE`, detailed progress messages will be printed to the console. Defaults to `FALSE`.
#' @param print_previous_ll_values If we want to print the previous values of the log-likelihood function. This can
#' be useful for controlling that the optimization procedure is proceeding in a monotone way and it does not
#' oscillate.
#' This argument is composed of two elements: TRUE/FALSE if we want or not to print the previous values and how many values we
#' want to print on the console. Default is (TRUE, 3), so that only the previous 3 values of the log-likelihood are printed.
#'
#' @return S3 object of class 'AdPaik', composed of several elements. See Details.
#' 
#'
#' @details The output of the model call 'AdPaikModel(...)' is a S3 object of class 'AdPaik', composed of the following quantities:
#' - formula: formula object provided in input by the user and specifying the relationship between the time-to-event, the covariates of
#' the dataset (regressors) and the cluster variable.
#' - dataset: matrix of the dataset containing the regressors and the dummy variables of the categorical covariates.
#' - Regressors: categorical vector of length R, with the name of the regressors.
#' They could be different from the original covariates of the dataset in case of categorical covariates.
#' Indeed, each categorical covariate with n levels needs to be transformed into (n-1) dummy variables and, therefore, (n-1) new regressors.
#' - NRegressors: number of regressors (R)
#' - ClusterVariable: name of the variable with respect to which the individuals can be grouped.
#' - NClusters: number of clusters/centres (also indicated with N).
#' - ClusterCodes: vector of length equal to the number of clusters, containing the codes of the clusters.
#' - TimeDomain: vector of the time-domain partition.
#' - NIntervals: number of intervals of the time-domain, also called with L. 
#' - NObservations: number of observations of the dataset.
#' - NParameters: number of parameters of the model. It can be computed as: \eqn{n_p = 2L + R + 2}.
#' - ParametersCategories: Numerical vector of length 5, containing the numerosity of each parameter category.
#' - ParametersRange: S3 object of class 'ParametersRange', containing ParametersRangeMin and ParametersRangeMax, two numerical vectors of length \eqn{n_p}, giving the minimum and the maximum range of each parameter, respectively.
#' - Loglikelihood: value of the maximized log-likelihood function, at the optimal estimated parameters.
#' - AIC: 'Akaike Information Criterion': it can be computed as \eqn{AIC = 2n_p - 2ll_{optimal}}.
#' It quantifies the loss of information related to the model fitting and output.
#' The smaller, the less the loss of information and the better the model accuracy.
#' - Status: Logical value. TRUE if the model reaches convergence, FALSE otherwise.
#' - NRun: Number of runs necessary to reach convergence. If the model does not reach convergence, such number is equal to the maximum number
#' of imposed runs.
#' - OptimalParameters: numerical vector of length \eqn{n_p}, containing the optimal estimated parameters (the parameters
#' that maximize the log-likelihood function).
#' - StandardErrorParameters: numerical vector of length \eqn{n_p}, corresponding to the standard error of each estimated parameters.
#' - ParametersCI: S3 object of class 'ParametersCI', composed of two numerical vector of length equal to \eqn{n_p}: the left and right 95\% confidence
#' interval of each estimated parameter of given level.
#' - BaselineHazard: numerical vector of length equal to L, containing the baseline hazard step-function.
#' - FrailtyDispersion:  S3 object of class 'FrailtyDispersion', containing two numerical vectors of length equal to L with the standard deviation and the variance of the frailty.
#' numerical vector of length equal to L (i.e. number of intervals of the time-domain), reporting the standard deviation
#' of the frailty.
#' - PosteriorFrailtyEstimates: S3 object of class 'PFE.AdPaik'. The object of class 'PFE.AdPaik' contains the Posterior Frailty Estimates computed with the procedure indicated in the reference paper and it is composed of three elements:
#'   - 'alpha': posterior frailty estimates for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
#'   - 'eps': posterior frailty estimates for \eqn{\epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'   - 'Z': posterior frailty estimates for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - PosteriorFrailtyVariance: S3 object of class 'PFV.AdPaik'. The object of class 'PFV.AdPaik' contains the Posterior Frailty Variances computed as indicated in the reference papaer and it is  composed of three elements:
#'   - 'alphaVar': posterior frailty variance for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
#'   - 'epsVar': posterior frailty variance for \eqn{\epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'   - 'ZVar': posterior frailty variance for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
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
#'\donttest{
#' # Call the main model
#' result <- AdPaikModel(formula, data_dropout, time_axis,
#'                       categories_range_min, categories_range_max, TRUE)
#' }

AdPaikModel <- function(formula, data, time_axis,
                        categories_range_min, categories_range_max,
                        n_extrarun = 60, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                        verbose = FALSE, print_previous_ll_values = c(TRUE, 3)){
  level = 0.95
  if (verbose) message("Adapted Paik et al.'s Model:")
  
  # Check all input variables are provided
  if(missing(categories_range_max))
    stop("At least one input variable is missing, with no default.")
  
  # Check time_axis vector
  check.time_axis(time_axis)
  
  # Check elements of the dataset
  check.dataset(data)
  
  # Check variables of formula are contained in dataset
  check.formula_terms(formula, data)
  
  # Extract elements of the formula object
  formula_variables <- all.vars(formula)
  
  # Extract position of cluster
  special <- c("cluster")
  terms_object <- terms(formula, specials = special, data = data)
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
  n_params <- 2 * n_intervals + n_regressors + 2
  n_run <- n_params + n_extrarun
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, 1, 1, n_intervals)
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
  
  # Change the warnings set to ignore warnings in the optimization phase
  # old_warnings <- getOption("warn")
  # suppressWarnings()
  
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
      result_optimize <- suppressWarnings(
        optimize(ll_AdPaik_1D,
                 c(params_range_min[index_to_vary], params_range_max[index_to_vary]),
                 maximum = TRUE, tol = tol_optimize,
                 index_to_vary, params, dataset, centre,
                 time_axis, dropout_matrix, e_matrix)
        )
      
      params[index_to_vary] <- result_optimize$maximum
    }
    
    global_optimal_params[r,] <- params
    global_optimal_loglikelihood_run <- ll_AdPaik_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix)
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
  
  # Set the warnings to the original value
  # options('warn' = old_warnings)
  
  # Extract best solution with maximum log-likelihood
  optimal_params <- global_optimal_params[optimal_run,]
  optimal_loglikelihood <- global_optimal_loglikelihood[optimal_run]
  
  # Compute the standard error from the Hessian matrix
  if (verbose) message(paste("Compute parameters standard error"))
  params_se <- params_se(optimal_params, params_range_min, params_range_max,
                         dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd)
  
  # Compute parameters confidence interval
  if (verbose) message(paste("Compute parameters confidence interval"))
  params_CI <- params_CI(optimal_params, params_se, level)
  
  # Compute baseline hazard step-function
  if (verbose) message(paste("Compute baseline hazard step function"))
  bas_hazard <- bas_hazard_internal(optimal_params, time_axis)
  
  # Compute frailty standard deviation
  if (verbose) message(paste("Compute frailty standard deviation"))
  frailty_dispersion <- frailty_Sd_internal(optimal_params, time_axis, n_regressors,
                                            categories_range_min, categories_range_max, TRUE)
  
  # Compute posterior frailty estimates
  if (verbose) message(paste("Compute posterior frailty estimates"))
  post_frailty_estimates <- post_frailty_internal(optimal_params, dataset, time_to_event, centre, time_axis)
  post_frailty_est <- post_frailty_estimates$PostFrailtyEst
  post_frailty_var <- post_frailty_estimates$PostFrailtyVar
  
  # Compute posterior frailty estimates confidence interval
  # if (verbose) message(paste("Compute posterior frailty estimates confidence interval"))
  # post_frailty_CI <- post_frailty_CI_internal(post_frailty_est, post_frailty_var, n_centres, n_intervals, level)
  
  # Akaike Information Criterium
  AIC = 2 * n_params - 2 * optimal_loglikelihood
  
  # Object to return
  return_list <- list("formula" = formula,
                      "dataset" = data[, formula_variables, drop = FALSE],
                      "Regressors" = new_covariates,
                      "NRegressors" = n_regressors,
                      "ClusterVariable" = cluster_name,
                      "NClusters" = n_centres,
                      "ClusterCodes" = centre_codes,
                      "TimeDomain" = time_axis,
                      "NIntervals" = n_intervals,
                      "NObservations" = n_individuals,
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
                      "PosteriorFrailtyEstimates" = post_frailty_est,
                      "PosteriorFrailtyVariance" = post_frailty_var)
  class(return_list) <- "AdPaik"
  
  # Return list of results
  return (return_list)
}


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
}#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
#' @title
#' One-Dimensional Analysis of log-Likelihood Function
#'
#' @description
#' Function for studying the log-likelihood function from the point of view of a single parameter and, therefore,
#' in a single direction.
#' It performs both the optimization of the log-likelihood with respect to this parameter and the
#' evaluation of the log-likelihood in several samples of the same parameter, while the other parameters can assume a constant assigned value or
#' can vary in their range.
#'
#' @details
#' The one-dimensional analysis of the log-likelihood function can be performed in two ways, with two different aims and results:
#' - Keeping fixed the other parameters (all the parameters in the vector, except for the one under consideration) to their optimal
#' value (flag_optimal_params = TRUE), determined through the multi-dimensional optimization. In this way, the optimized value of the parameter
#' coincides with the one get with the general and global approach and, therefore, there is no need to repeat this procedure several times (n_iter = 1).
#' However, this approach is really useful if we want to check the trend the log-likelihood function and to observe if it increases, decreases or
#' is constant.
#' - Letting the other parameters vary in their range (flag_optimal_params = FALSE). The optimized parameter value will always assume a different
#' value, because it depends on the value of the other parameters, and it is suggested to repeat the procedure several times (n_iter \eqn{\geq} 5),
#' so that it is possible to identify a precise existence region for such parameter
#'
#' @param formula Formula object indicating the response variable, the covariates and the cluster variable.
#' @param data Dataset in which the variables of the formula object are located.
#' @param time_axis Partitioned time-domain.
#' @param index_param_to_vary Index of the parameter, in the parameter vector, with respect to which the log-likelihood function
#' is maximized in a one-dimensional way. The index s provided to identify the parameter under consideration inside the vector, avoiding
#' providing its name or value.
#' @param flag_optimal_params Are the other parameters extracted from the optimal vector of parameters? If so, the flag should be equal to TRUE.
#' Otherwise, the flag is equal to FALSE.
#' @param optimal_params Vector of optimal parameters, determined through an entire multi-dimensional maximization
#' of the log-likelihood function. The default value (NULL) indicates that no vector is provided
#' and the parameters are randomly extracted in their range.
#' @param categories_range_min Vector containing the minimum value assumed by each parameter category.
#' @param categories_range_max Vector containing the maximum value assumed by each parameter category.
#' @param n_iter Number of times the one-dimensional analysis with respect to the indicated parameter must be executed.
#' Default value is 5. See details for more information.
#' @param tol_optimize Tolerance used in the optimize R function for the one-dimensional optimization
#' of the log-likelihood function.
#' @param flag_plot Logical value for plotting the trend of the log-likelihood function with respect to the parameter under consideration.
#' A plot for each iteration (n_iter) is reported. Defaults to FALSE.
#' Be careful that if the optimal parameters are provided, then the trend may be always the same and therefore it may be sufficient to
#' set n_iter = 1. On the other hand, if optimal parameters are not provided, then it is recommended to impose a higher n_iter.
#' @param n_points Number of internal points in which the log-likelihood function must be evaluated, to
#' plot it.
#' @param cex Dimension of the points in the plot.
#' @param cex_max Dimension of the optimal point in the plot.
#' @param color_bg Color used in the plot for the points.
#' @param color_max_bg Color used for the optimal point in the plot.
#' @param pch Shape to be used for the points.
#'
#' @return If the flag for the plot has been activated, the function returns both the plot of the one-dimensional log-likelihood function and
#' a class S3 object. Otherwise, only a S3 object of class 'AdPaik_1D'.
#' This class object is composed of:
#' - numerical vector of length @n_iter containing the optimal estimated parameter.
#' - numerical vector of length @n_iter containing the associated one-dimensional optimized log-likelihood value
#'
#' @importFrom stats terms runif optimize
#' @importFrom graphics lines points legend
#' 
#' @export
#'
#' @examples
#' # Consider the 'Academic Dropout dataset'
#' data(data_dropout)
#' # Define the variables needed for the model execution
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' eps <- 1e-10
#' \donttest{
#' # Identify a parameter existence range
#' categories_range_min <- c(-8, -2, eps, eps, eps)
#' categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)
#' index_param_to_vary <- 1
#' analysis_1D_opt <- AdPaik_1D(formula, data_dropout,
#'                              time_axis, index_param_to_vary, 
#'                              flag_optimal_params = FALSE, 
#'                              optimal_params = NULL,
#'                              flag_plot = TRUE,
#'                              categories_range_min, categories_range_max, 
#'                              n_iter = 5)
#' 
#'
#' # or Study the log-likelihood behaviour
#' categories_range_min <- c(-8, -2, eps, eps, eps)
#' categories_range_max <- c(-eps, 0.4, 1 - eps, 1, 10)
#' index_param_to_vary <- 14
#' # Call the main model
#' result <- AdPaikModel(formula, data_dropout, time_axis,
#'                       categories_range_min, categories_range_max, TRUE)
#' analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
#'                              index_param_to_vary, flag_optimal_params = TRUE, 
#'                              flag_plot = TRUE, optimal_params = result$OptimalParameters,
#'                              categories_range_min, categories_range_max, n_iter = 1)
#' }


AdPaik_1D <- function(formula, data, time_axis,
                      index_param_to_vary, flag_optimal_params = FALSE, 
                      optimal_params = NULL,
                      categories_range_min, categories_range_max,
                      n_iter = 5, tol_optimize = 1e-6,
                      flag_plot = FALSE, n_points = 150,
                      cex = 0.7, cex_max = 0.8, 
                      color_bg = "black", color_max_bg = "red",
                      pch = 21){
  
  # Check all input variables are provided
  if(missing(categories_range_max))
    stop("At least one input variable is missing, with no default.")
  
  # Check time_axis vector
  check.time_axis(time_axis)
  
  # Check elements of the dataset
  check.dataset(data)
  
  # Check variables of formula are contained in dataset
  check.formula_terms(formula, data)
  
  # Check either the optimal parameters are provided or they must be simulated
  check.flag_optimal_params(optimal_params, flag_optimal_params)
  
  # Extract elements of the formula object
  formula_variables <- all.vars(formula)
  
  # Extract position of cluster
  special <- c("cluster")
  terms_object <- terms(formula, specials = special, data = data)
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
      dummy_extracted <- extract_dummy_variables(data[,covariates[j]], 
                                                 covariates[j])
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
  n_params <- 2 * n_intervals + n_regressors + 2
  
  # Check existence provided index
  check.index(index_param_to_vary, n_params)
  
  # Define vector of categories for Adapted Paik et al.'s Model
  params_categories <- c(n_intervals, n_regressors, 1, 1, n_intervals)
  n_categories <- length(params_categories)
  
  # Check the correctness of provided range categories
  check.categories_params(n_categories, categories_range_min, categories_range_max)
  
  # Generate extended vector of parameters ranges
  params_range_min <- params_range_max <- c()
  for(c in 1: n_categories){
    n_params_in_c <- params_categories[c]
    params_range_min <- c(params_range_min, 
                          rep(categories_range_min[c], n_params_in_c))
    params_range_max <- c(params_range_max, 
                          rep(categories_range_max[c], n_params_in_c))
  }
  
  # Check that provided optimal parameters are contained in their min, max range
  if(flag_optimal_params)
    check.range_params(optimal_params, params_range_min, params_range_max)
  
  # Build the matrices e_{ijk} and d_{ijk}
  e_matrix <- matrix(rep(0, n_intervals * n_individuals), 
                     n_individuals, n_intervals)
  dropout_matrix <- matrix(rep(0, n_intervals * n_individuals), 
                           n_individuals, n_intervals)
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
      params[index_param_to_vary] <- runif(1, params_range_min[index_param_to_vary], 
                                           params_range_max[index_param_to_vary])
    }
    
    # Optimize the log-likelihood wrt the indicated parameter
    result_optimize <- suppressWarnings(
      optimize(ll_AdPaik_1D,
               c(params_range_min[index_param_to_vary], 
                 params_range_max[index_param_to_vary]),
               maximum = TRUE, tol = tol_optimize,
               index_param_to_vary, params, dataset, centre,
               time_axis, dropout_matrix, e_matrix)
    )
    
    param_optimal[iter] <- result_optimize$maximum
    ll_optimized[iter] <- result_optimize$objective
    
    if(flag_plot){
      param_1D <- param_optimal[iter]
      ll_1D <- ll_optimized[iter]
      plot_ll_1D(param_1D, index_param_to_vary, ll_1D, params,
                 params_range_min[index_param_to_vary], 
                 params_range_max[index_param_to_vary],
                 dataset, centre, time_axis, dropout_matrix, e_matrix,
                 n_points, cex, cex_max, color_bg, color_max_bg, pch)
    }
  }
  
  return_list <- list("EstimatedParameter" = param_optimal,
                      "OptimizedLoglikelihood" = ll_optimized)
  
  return (return_list)
}
