#' @title
#' Adapted Paik et Al.'s Model: Time-Dependent Shared Frailty Cox Model
#'
#' @description
#' Function for applying the 'Adapted Paik et al.'s Model', an innovative Cox Model with time-dependent frailty,
#' shared by individuals belonging to the same group/cluster.
#'
#' To generate time-dependence, the temporal domain must be divided into a certain number
#' of intervals. For more information about the time-domain, its relationship with the follow-up and its internal
#' subdivision refer to Details.
#'
#' The model log-likelihood function depends on a certain number of parameters and it is maximized with respect to all of them,
#' using a reinterpretation of the 'Powell's method in multidimension', that is a multi-dimensional optimization method based on
#' repeated one-dimensional optimization of the log-likelihood function (with respect to one parameter at the time).
#' In this context, the one-dimensional optimization is performed through the 'optimize' R function.
#' For more information about the unknown model parameters, their type and numerosity refer to Details.
#'
#' Several quantities are estimated at the end of the optimization phase:
#' - parameters, their standard error and confidence interval;
#' - baseline hazard;
#' - frailty dispersion (standard deviation and variance);
#' - posterior frailty estimates, their variance and confidence interval;
#' - Akaike Information Criterion (AIC).
#'
#' @details
#' Two observation needs to made about the time-domain:
#' - The time domain may coincide with the follow-up or it may be contained in it.
#' Indeed, the left boundary can be greater than the beginning of the follow-up and, for instance, it can coincide
#' with the time-instants in which the events begin to happen; conversely, the right boundary of the two must be the same.
#' - The partition of the time domain into intervals can be made according to two selected criteria:
#' (1) using an already existent partition of the follow-up
#' (2) using the shape of the baseline hazard function as reference: divide the time-domain according to regions in
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
#' The category vector is structured as follows: (baseline log-hazard, regressors, \eqn{\mu_1}, \eqn{\nu}, \eqn{\gamma_k}) with dimension
#' (n_intervals, n_regressors, 1, 1, n_intervals).
#'
#'
#' @param formula Formula object having on the left hand side the @time_to_event variable, that is the time-instant in which
#' the individual failed. On the right hand side, it has the regressors and the cluster variable.
#' @param data Dataset in which all variables of the formula object must be found and contained.
#' This dataset can also contain other variables, but they will not be considered.
#' It can be either a dataframe or a matrix, but in both cases the name of each column must correspond to the name of
#' the formula variables. It is not necessary to attach it (in case of a data.frame)
#' @param time_axis Temporal domain
#' @param categories_range_min Vector containing the minimum value assumable by each parameter category.
#' @param categories_range_max Vector containing the maximum value assumable by each parameter category.
#' @param flag_fullsd Logical. If TRUE, the full frailty standard deviation is computed, otherwise the partial one that keeps into
#' account only the time-dependent component. Defaults to `TRUE`.
#' @param n_extrarun Total number of runs (iterations) are obtained summing to the number of parameters and n_extrarun.
#' @param tol_ll Tolerance on the log-likelihood value.
#' @param tol_optimize Internal tolerance for the one-dimensional optimization through 'optimize' R function.
#' @param h_dd Discretization step used for the numerical approximation of the second derivative of the log-likelihood function.
#' @param print_previous_ll_values If we want to print the previous values of the log-likelihood function. This can
#' be useful for controlling that the optimization procedure is proceeding in a monotone way and it does not
#' oscillate.
#' This argument is composed of two elements: TRUE/FALSE if we want or not to print the previous values and how many values we
#' want to print on the console. Default is (TRUE, 3), so that only the previous 3 values of the log-likelihood are printed.
#' @param level A numeric value representing the confidence level for the optimal parameters (default is 0.95 for 95% confidence).
#' @param verbose Logical. If `TRUE`, detailed progress messages will be printed to the console. Defaults to `FALSE`.
#'
#' @return S3 object of class 'AdPaik', composed of several elements. See Details.
#'
#' @details The output of the model call 'AdPaikModel(...)' is a S3 object of class 'AdPaik', composed of the following quantities:
#' - formula: formula object provided in input by the user and specifying the relationship between the time-to-event, the covariates of
#' the dataset (regressors) and the cluster variable.
#' - Regressors: categorical vector of length R, with the name of the regressors.
#' They could be different from the original covariates of the dataset in case of categorical covariates.
#' Indeed, each categorical covariate with n levels needs to be transformed into (n-1) dummy variables and, therefore, (n-1) new regressors.
#' - NRegressors: number of regressors (R)
#' - ClusterVariable: name of the variable with respect to which the individuals can be grouped.
#' - NClusters: number of clusters/centres (also indicated with N).
#' - ClusterCodes: vector of length N, containing the codes of the clusters.
#' - NIntervals: number of intervals of the time-domain, also called with L. 
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
#' - PosteriorFrailtyEstimates: S3 object of class 'PFE.AdPaik'. See details.
#' - PosteriorFrailtyVariance: S3 object of class 'PFV.AdPaik'. See details.
#' - PosteriorFrailtyCI: S3 object of class 'PFCI.AdPaik'. See details.
#'
#' @details
#' The object of class 'PFE.AdPaik' contains the Posterior Frailty Estimates computed with the procedure indicated in the reference paper and
#' it is composed of three elements:
#' - 'alpha': posterior frailty estimates for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
#' - 'eps': posterior frailty estimates for \eqn{\epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - 'Z': posterior frailty estimates for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'
#' @details
#' The object of class 'PFV.AdPaik' contains the Posterior Frailty Variances computed as indicated in the reference papaer and it
#' is  composed of three elements:
#' - 'alphaVar': posterior frailty variance for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
#' - 'epsVar': posterior frailty variance for \eqn{\epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - 'ZVar': posterior frailty variance for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'
#' @details
#' The object of class 'PFCI.AdPaik' contains the Posterior Frailty Confidence Interval and it is composed of two elements:
#' - left confidence interval for the estimated \eqn{\hat{Z}_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - right confidence interval for the estimated \eqn{\hat{Z}_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'
#' @source
#' ...
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
                        flag_fullsd = TRUE,
                        n_extrarun = 60, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                        print_previous_ll_values = c(TRUE, 3),
                        level = 0.95,
                        verbose = FALSE){
  
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
  bas_hazard <- bas_hazard(optimal_params, time_axis)
  
  # Compute frailty standard deviation
  if (verbose) message(paste("Compute frailty standard deviation"))
  frailty_dispersion <- frailty_Sd.AdPaik(optimal_params, time_axis, n_regressors,
                                          categories_range_min, categories_range_max, TRUE)
  
  # Compute posterior frailty estimates
  if (verbose) message(paste("Compute posterior frailty estimates"))
  post_frailty_estimates <- post_frailty.AdPaik(optimal_params, dataset, time_to_event, centre, time_axis)
  post_frailty_est <- post_frailty_estimates$PostFrailtyEst
  post_frailty_var <- post_frailty_estimates$PostFrailtyVar
  
  # Compute posterior frailty estimates confidence interval
  if (verbose) message(paste("Compute posterior frailty estimates confidence interval"))
  post_frailty_CI <- post_frailty_CI.AdPaik(post_frailty_est, post_frailty_var, n_centres, n_intervals, level)
  
  # Akaike Information Criterium
  AIC = 2 * n_params - 2 * optimal_loglikelihood
  
  # Object to return
  return_list <- list("formula" = formula,
                      "dataset" = data[, formula_vars, drop = FALSE],
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
                      "PosteriorFrailtyVariance" = post_frailty_var,
                      "PosteriorFrailtyCI" = post_frailty_CI)
  class(return_list) <- "AdPaik"
  
  # Return list of results
  return (return_list)
}

