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
  class(return_list) <- "AdPaik_1D"
  
  return (return_list)
}