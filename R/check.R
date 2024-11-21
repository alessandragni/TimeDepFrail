#' @title
#' Check correctness of time domain subdivision
#'
#' @description
#' The function controls that the time domain is a vector and it has at least 2 elements and that all of them are not negative.
#' Moreover, it checks that all its elements are non-negative and properly ordered, in an ascending way.
#'
#'
#' @param time_axis Numerical vector of temporal domain subdivision.
#'
#' @return An error is returned if any condition is not satisfied.

check.time_axis <- function(time_axis){
  
  # Check structure of time_axis
  if(! is.vector(time_axis))
    stop("'time_axis' input variable must be a vector.")
  
  # Check all elements of the time_axis vector are non-negative
  # Check length is greater tha 1
  length_time_axis <- length(time_axis)
  
  if(length_time_axis <= 1)
    stop("'time_axis' input variable must have at least 2 elements.")
  
  for(k in 1:length_time_axis){
    if(time_axis[k] < 0)
      stop("Negative element in 'time_axis' input variable.")
  }
  
  # Check the time_axis vector is properly ordered
  first_time <- time_axis[1]
  for(k in 2:length_time_axis){
    if(first_time >= time_axis[k])
      stop("Not ordered 'time_axis' input variable.")
    first_time <- time_axis[k]
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check structure of Posterior Frailty Estimates
#'
#' @description
#' The function controls that the structure of the 'Posterior Frailty Estimates' coincides with the theoretical one.
#'
#' @param post_frailty_est Posterior frailty estimates S3 object of class 'PFE.AdPaik', composed of three elements:
#' - 'alpha': posterior frailty estimates for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
#' - 'eps': posterior frailty estimates for \eqn{\epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals).
#' - 'Z': posterior frailty estimates for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals)
#' @param n_intervals Number of intervals of the time-domain
#' @param n_centres Number of centres/clusters.
#'
#' @return An error if any condition is not satisfied.

check.structure_post_frailty_est <- function(post_frailty_est, n_intervals, n_centres){
  # Check class
  if (!inherits(post_frailty_est, "PFE.AdPaik"))
    stop("First argument is not of class 'PFE.AdPaik'.")
  
  if((! is.list(post_frailty_est)) || (length(post_frailty_est) != 3))
    stop("Wrong structure for first input variable 'post_frailty_est'.")
  
  if(names(post_frailty_est)[1] != "alpha")
    stop("Wrong name for the first element of 'post_frailty_est' variable.")
  if(! is.vector(post_frailty_est$alpha))
    stop("Wrong structure for the first element of 'post_frailty_est' variable.")
  if(length(post_frailty_est$alpha) != n_centres)
    stop("Different values for the number of centres and number of time-independent posterior frailty estimates (alpha).")
  
  if(names(post_frailty_est)[2] != "eps")
    stop("Wrong name for the second element of 'post_frailty_est' variable.")
  if(! is.matrix(post_frailty_est$eps))
    stop("Wrong structure for the second element of 'post_frailty_est' variable.")
  if((dim(post_frailty_est$eps)[1]) != n_centres)
    stop("Different values for the number of centres and the time-varying posterior estimates (eps).")
  if((dim(post_frailty_est$eps)[2]) != n_intervals)
    stop("Different values for the number of intervals and the time-varying posterior estimates (eps).")
  
  if(names(post_frailty_est)[3] != "Z")
    stop("Wrong name for the third element of 'post_frailty_est' variable.")
  if(! is.matrix(post_frailty_est$Z))
    stop("Wrong structure for the third element of 'post_frailty_est' variable.")
  if((dim(post_frailty_est$Z)[1]) != n_centres)
    stop("Different values for the number of centres and the overall time-varying posterior estimates.")
  if((dim(post_frailty_est$Z)[2]) != n_intervals)
    stop("Different values for the number of intervals and the overall time-varying posterior estimates.")
}
#-------------------------------------------------------------------------------
#' @title
#' Check structure of Posterior Frailty Variances
#'
#' @description
#' The function controls that the structure of the 'Posterior Frailty Variances' coincides with the theoretical one.
#'
#' @param post_frailty_var Posterior frailty variances S3 object of class 'PFV.AdPaik', composed of three elements:
#' - 'alphaVar': posterior frailty variance for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of centres.
#' - 'epsVar': posterior frailty variance for \eqn{\epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals).
#' - 'ZVar': posterior frailty variance for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. It is a matrix of dimension (n_centres, n_intervals)
#' @param n_intervals Number of intervals of the time-domain
#' @param n_centres Number of centres/clusters.
#'
#' @return An error if any condition is not satisfied.
#' 
check.structure_post_frailty_var <- function(post_frailty_var, n_intervals, n_centres){
  # Check class
  if(!inherits(post_frailty_var, "PFV.AdPaik"))
    stop("First argument is not of class 'PFV.AdPaik'.")
  
  if((! is.list(post_frailty_var)) || (length(post_frailty_var) != 3))
    stop("Wrong structure for first input variable (post_frailty_var).")
  
  if(names(post_frailty_var)[1] != "alphaVar")
    stop("Wrong name for the first element of 'post_frailty_var' variable.")
  if(! is.vector(post_frailty_var$alphaVar))
    stop("Wrong structure for the first element of 'post_frailty_var' variable.")
  if(length(post_frailty_var$alpha) != n_centres)
    stop("Different values for the number of centres and number of time-independent posterior frailty variances (alphaVar).")
  
  if(names(post_frailty_var)[2] != "epsVar")
    stop("Wrong name for the second element of 'post_frailty_var' variable.")
  if(! is.matrix(post_frailty_var$epsVar))
    stop("Wrong structure for the second element of 'post_frailty_var' variable.")
  if((dim(post_frailty_var$epsVar)[1]) != n_centres)
    stop("Different values for the number of centres and the time-varying posterior variances (epsVar).")
  if((dim(post_frailty_var$epsVar)[2]) != n_intervals)
    stop("Different values for the number of intervals and the time-varying posterior variances (epsVar).")
  
  if(names(post_frailty_var)[3] != "ZVar")
    stop("Wrong name for the third element of 'post_frailty_var' variable.")
  if(! is.matrix(post_frailty_var$ZVar))
    stop("Wrong structure for the third element of 'post_frailty_var' variable.")
  if((dim(post_frailty_var$Z)[1]) != n_centres)
    stop("Different values for the number of centres and the overall time-varying posterior variances.")
  if((dim(post_frailty_var$Z)[2]) != n_intervals)
    stop("Different values for the number of intervals and the overall time-varying posterior variances.")
}
#-------------------------------------------------------------------------------
#' @title
#' Check non-negativeness of the posterior frailty estimates
#'
#' @description
#' The function controls that all posterior frailty estimates are non-negative.
#' Indeed, by construction the realizations of a gamma distribution are non negative.
#'
#' @param post_frailty_est An S3 class object containing the posterior frailty estimates: \eqn{\hat{\alpha}_j, \hat{\epsilon}_{jk},
#'  \hat{Z}_{jk}, \forall j,k}
#' @param n_centres Number of groups/clusters.
#' @param n_intervals Number of intervals of the time domain
#'
#' @return An error if any condition is not satisfied.

check.value_post_frailty <- function(post_frailty_est, n_centres, n_intervals){
  # Check that no estimates are negative
  for (i in 1:n_centres){
    if(post_frailty_est$alpha[i] < 0)
      stop("Posterior frailty estimate for 'alpha' must be non-negative.")
  }
  
  for(i in 1:n_centres){
    for(k in 1:n_intervals){
      if(post_frailty_est$eps[i,k] < 0)
        stop("Posterior frailty estimate for 'eps' must be non negative.")
      if(post_frailty_est$Z[i,k] < 0)
        stop("Posterior frailty estimates for 'Z' must be non negative.")
    }
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check numerosity of posterior frailty estimates
#'
#' @description
#' The function controls that a time-dependent posterior frailty estimate is computed for each centre
#'
#' @param post_frailty_est An S3 class object containing the posterior frailty estimates \eqn{\hat{\alpha}_j, \hat{\epsilon}_{jk},
#'  \hat{Z}_{jk}, \forall j,k}
#' @param centre_codes Numerical vector of length equal to the number of distinct centres/clusters in the study
#'
#' @return An error if any condition is not satisfied.

check.post_frailty_centre <- function(post_frailty_est, centre_codes){
  # Extract number of centres
  n_centres <- length(centre_codes)
  
  # Check that the dimension of the post frailty estimate coincide with the number of centre
  if(length(post_frailty_est$alpha) != n_centres)
    stop("Different values for the number of centres and the number of posterior frailty estimates (alpha).")
  
  if(dim(post_frailty_est$eps)[1] != n_centres)
    stop("Different values for the number of centres and the number of posterior frailty estimates (eps).")
}
#-------------------------------------------------------------------------------
#' @title
#' Check structure of Posterior Frailty Confidence Interval
#'
#' @description
#' The function controls that the structure of the 'Posterior Frailty Confidence Interval' coincides with the theoretical one.
#'
#' @param post_frailty_CI Posterior frailty estimates S3 object of class 'PFCI.AdPaik', composed of two elements:
#' - left confidence interval for the estimated \eqn{\hat{Z}_{jk}, \forall j,k}
#' - right confidence interval for the estimated \eqn{\hat{Z}_{jk}, \forall j,k}
#' @param n_intervals Number of intervals of the time-domain
#' @param n_centres Number of centres/clusters.
#'
#' @return An error if any condition is not satisfied.
check.structure_post_frailty_CI <- function(post_frailty_CI, n_intervals, n_centres){
  # Check class
  if(!inherits(post_frailty_CI, "PFCI.AdPaik"))
    stop("First argument is not of class 'PFCI.AdPaik'.")
  
  # Check structure
  if((! is.list(post_frailty_CI)) || (length(post_frailty_CI) != 2))
    stop("Wrong structure for first input variable (post_frailty_CI).")
  
  # Check internal structures
  if(names(post_frailty_CI)[1] != "PostFrailtyCI_left")
    stop("Wrong name for the first element of 'post_frailty_CI' variable.")
  if(! is.matrix(post_frailty_CI$PostFrailtyCI_left))
    stop("Wrong structure for the first element of 'post_frailty_CI' variable.")
  if((dim(post_frailty_CI$PostFrailtyCI_left)[1]) != n_centres)
    stop("Different values for the number of centres and the left confidence interval.")
  if((dim(post_frailty_CI$PostFrailtyCI_left)[2]) != n_intervals)
    stop("Different values for the number of intervals and the left confidence interval.")
  
  if(names(post_frailty_CI)[2] != "PostFrailtyCI_right")
    stop("Wrong name for the second element of 'post_frailty_CI' variable.")
  if(! is.matrix(post_frailty_CI$PostFrailtyCI_right))
    stop("Wrong structure for the second element of 'post_frailty_CI' variable.")
  if((dim(post_frailty_CI$PostFrailtyCI_right)[1]) != n_centres)
    stop("Different values for the number of centres and the right confidence interval.")
  if((dim(post_frailty_CI$PostFrailtyCI_right)[2]) != n_intervals)
    stop("Different values for the number of intervals and the right confidence interval.")
}
#-------------------------------------------------------------------------------
#' @title
#' Check correctness for the cluster variable
#'
#' @description
#' The function controls that the provided cluster variable is a vector, with at least two levels. Indeed, it is not possible to
#' apply the Time-Dependent Shared Frailty Cox Model with no clusters.
#'
#' @param centre Numerical vector of length equal to the number of individuals in the study, containing the individual grouo/cluster
#' membership.
#'
#' @return An error if any condition is not satisfied

check.centre <- function(centre){
  # Check that it is a vector
  if(! is.vector(centre))
    stop("Wrong structure for 'cluster' variable.")
  
  # Extract levels from input variable and compute number of centres
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)
  
  # At least two centre
  if(n_centres == 0)
    stop("Null dimension for the 'cluster' variable.")
  else if(n_centres == 1)
    stop("Only one cluster in the 'cluster' variable.")
}

#-------------------------------------------------------------------------------
#' @title
#' Check correctness of plot variables pch and color
#'
#' @description
#' The function controls that the input variables 'pch_type' and 'color_bg' have the correct structure,
#' they have the same dimension of the number of clusters in the dataset and they have meaningful elements.
#'
#' These variables are used for the plot of the posterio frailty estimates: the estimates for each faculty are plotted through a symbol, having
#' color and shape indicated by the variables (for the k-th faculty, consider the k-th element of both vectors).
#'
#' @param centre_codes Numerical vector of length equal to the number of centres/clusters in the dataset and containing the distinct centres/clusters.
#' They correspond to the levels of the numerical vector of individual group membership.
#' @param pch_type Numerical vector of length equal to the number of centres and containing the point shape for each faculty.
#' @param color_bg Numerical vector of length equal to the number of centres and containing the color of the point for each faculty.
#'
#' @return An error if any condition is not satisfied.

check.pchtype_colorbg <- function(centre_codes, pch_type, color_bg){
  # Check structure
  if(! is.vector(pch_type))
    stop("Wrong structure for pch_type variable.")
  if(! is.vector(color_bg))
    stop("Wrong structure for color_bg variable.")
  
  # Compute number of centres and check it coincide with the dimension of the two provided vectors
  n_centres <- length(centre_codes)
  if(length(pch_type) != n_centres)
    stop("Wrong number of element for the pch_type variable. It does not correspond to the number of different centres.")
  
  if(length(color_bg) != n_centres)
    stop("Wrong number of element for the color_bg variable. It does not correspond to the number of different centres.")
  
  # Check that elements of the two vectors are properly assigned
  for(i in 1:n_centres){
    if(! is.character(color_bg[i]))
      stop("Element of color_bg vector variable is not a string.")
    
    if((pch_type[i] < 0) || (pch_type[i] > 25))
      stop("Element of pch_type vector variable outside the range 0-25.")
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check correctness of legend position
#'
#' @description
#' The function controls that the provided position of the legend is correct. It can be either
#' a vector of length 2, giving the x and y coordinates, or a string, giving the exact position among different possibilities.
#'
#' @param pos_legend Either a numerical vector of length 2, with the x and y coordinates, or a string with the exact position.
#'
#' @return An error if any condition is not satisfied.

check.poslegend <- function(pos_legend){
  
  # Check it is either a vector or a string
  if((!is.vector(pos_legend)) & (! is.character(pos_legend)))
    stop("Wrong pos_legend variable.")
  
  # In case it is a vector, it must have two elements
  if(! is.character(pos_legend)){
    if(length(pos_legend) != 2)
      stop("pos_legend vector variable must have two elements (x, y).")
  }
  
  # In case it is a string, look inside the admissible values
  possible_pos_legend <- c("bottomright", "bottom", "bottomleft", "left", "topleft", "top",
                           "topright", "right", "center")
  n_possible_pos_legend <- length(possible_pos_legend)
  counter <- 0
  if(is.character(pos_legend)){
    for(i in 1:n_possible_pos_legend){
      if(pos_legend == possible_pos_legend[i])
        counter <- counter + 1
    }
    if(counter != 1)
      stop("Wrong pos_legend provided")
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check correctness of frailty standard deviation
#'
#' @description
#' The function controls that the frailty standard deviation vector has a length equal to the number of inyervals of the time domain and
#' that its elements are non-negative.
#'
#' @param frailty_dispersion Frailty dispersion
#' @param n_intervals Number of intervals of the time-domain
#'
#' @return An error if any condition is not satisfied.

check.frailty_dispersion <- function(frailty_dispersion, n_intervals){
  
  # Check structure
  if(!inherits(frailty_dispersion, 'FrailtyDispersion'))
    stop("Wrong class of 'frailty_dispersion'.")
  
  if(length(frailty_dispersion) != 2)
    stop("Wrong number of elements in 'frailty_dispersion' class.")
  
  # Check input structure of sd
  if(! is.vector(frailty_dispersion$FrailtyStandardDeviation))
    stop("Frailty standard deviation is not a vector.")
  if(! is.vector(frailty_dispersion$FrailtyVariance))
    stop("Frailty Variance is not a vector.")
  
  # Check length of the sd vector
  if(n_intervals != length(frailty_dispersion$FrailtyStandardDeviation))
    stop("Length of standard deviation vector different from number of intervals.")
  if(n_intervals != length(frailty_dispersion$FrailtyVariance))
    stop("Length of standard deviation vector different from number of intervals.")
}
#-------------------------------------------------------------------------------
#' @title
#' Check positiviness of the frailty standard deviation
#' 
#' @description
#' The method controls that the frailty standard deviation vector has non-negative elements
#' 
#' @param sd Numerical vector of length equal to the number of intervals, containing the frailty standard deviation
#' @param n_intervals Number of intervals of the time-domain
#' 
#' @return An error if any condition is not satisfied.

check.pos_frailty_sd <- function(sd, n_intervals){
  # Check correctness of dimension of sd vector
  if(length(sd) != n_intervals)
    stop("Length of standard deviation vector different from number of intervals.")
  
  # Check positiviness of sd vector
  for(k in 1:n_intervals){
    if(sd[k] < 0)
      stop("Negative value for the frailty standard deviation.")
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check correctness of formula terms
#'
#' @description
#' The function controls that the terms composing the formula object provided in input to the main model are correct.
#' They must include:
#' - response variable on the left hand side
#' - covariates (numerical or categorical) on the right hand side
#' - cluster variable (categorical) on the right hand side and specified by the function 'cluster()'
#'
#' Moreover, it controls that the covariates are contained in the dataset provided.
#'
#' @param formula Formula object specifying the relationship between the time-to-event, the covariates and the cluster variables.
#' @param data Dataset in which these variables can be found.
#'
#' @return An error if any condition is not satified.

check.formula_terms <- function(formula, data){
  # Extract terms from the formula object
  special <- c("cluster")
  terms_object <- terms(formula, specials = special, data = data)
  
  # Check response
  response <- attr(terms_object, "response")
  if( response != 1)
    stop("No response (left hand side) in the formula object.")
  
  # Check covariates
  covariates <- attr(terms_object, "term.labels")
  n_covariates <- length(covariates)
  if(n_covariates == 0)
    stop("No covariates in the formula object.")
  
  # Check cluster
  cluster_index <- attr(terms_object, "specials")$cluster
  if(is.null(cluster_index))
    stop("No 'cluster' variable in the formula object.")
  check.centre(data[,cluster_index])
  
  # Extract number of covariates
  n_covariates <- n_covariates - 1
  
  # Define the number of variables of the formula object
  n_variables <- n_covariates + 2
  formula_variables <- all.vars(terms_object)
  if(length(formula_variables) != n_variables)
    stop("Wrong computed number of variables in the formula object.")
  
  # Extract the name and numerosity of the dataset columns
  ncol_data <- dim(data)[2]
  data_names <- colnames(data)
  if(n_variables > ncol_data)
    stop("Number of formula variables higher than number of dataset columns.")
  
  # Check formula variables contained in data
  for(i in 1:n_variables){
    counter <- 0
    for(j in 1:ncol_data){
      if(formula_variables[i] == data_names[j])
        counter <- 1
    }
    if(counter != 1){
      string_error <- paste("Variable '", formula_variables[i], "' not contained in the dataset.")
      stop(string_error)
    }
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check correctness of parameters categories
#'
#' @description
#' The function controls that the provided parameters categories have a length equal to the number of categories required by the model
#' parameters. For the current model, the number of categories is 5 because there are five blocks of unkown parameters (\eqn{\phi_k \forall k, \beta_r \forall r,
#' \mu_1, \nu, \gamma_k \forall k}).
#'
#' Moreover, it also controls that the minimum value of a parameter category is actually less than or eqaul to the maximum value for the same category.
#'
#' @param n_categories Number of categories expected by the model. For the current model, they are 5.
#' @param categories_range_min Numerical vector of length 5, containing the minimum ranges for the parameters beloning to those categories.
#' @param categories_range_max Numerical vector of length equal to 5, containing the maximum ranges for the parameters belonging to those categories.
#'
#' @return An error if the any condition is not satisfied.

check.categories_params <- function(n_categories, categories_range_min, categories_range_max){
  # Check dimension of categories vectors
  if((length(categories_range_min) != n_categories) || (length(categories_range_max) != n_categories))
    stop("Provided wrong length of either minimum or maximum categories range.")
  
  # Check the minimum value indicated is less than or equal to the maximum value
  for(p in 1:n_categories){
    if(categories_range_min[p] > categories_range_max[p])
      stop("Minimum value for a category > maximum value for the same category.")
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check correctness of input parameters
#'
#' @description
#' The function controls that the input parameter vector have a length equal to the theoretical one required by the model and that each parameter
#' properly belongs to its range.
#'
#' @param optimal_params Numerical vector of length equal to the number of model parameters. For the 'Adapted Paik et al.'s Model' it can be computed
#' as: \eqn{n_p = 2L + R + 2}, where \eqn{L} stands for the number of intervals of the time domain and \eqn{R} the number of regressors of the dataset.
#' @param params_range_min Numerical vector of length equal to the number of model parameters (\eqn{n_p}) and containing the minimum range for each
#' parameter.
#' @param params_range_max Numerical vector of length equal to the number of model parameters (\eqn{n_p}) and containing the maximum range for each
#' parameter.
#'
#' @return An error if any condition is not satisfied.

check.range_params <- function(optimal_params, params_range_min, params_range_max){
  # Extract numerosity of the parameters
  n_params <- length(optimal_params)
  
  # Check that actual numerosity coincides with the expected one
  if(length(params_range_min) != n_params)
    stop("Different number of elements in 'optimal_params' and 'params_range_min'.")
  if(length(params_range_max) != n_params)
    stop("Different number of elements in 'optimal_params' and 'params_range_max'.")
  
  # Check that each parameter belongs to its range
  for(p in 1:n_params){
    if((optimal_params[p] < params_range_min[p]) || (optimal_params[p] > params_range_max[p])){
      error_string <- paste("Parameter number ", p, "not in its range.")
      stop(error_string)
    }
  }
}
#-------------------------------------------------------------------------------
check.time_domain <- function(time_domain, flag_time_domain){
  if(flag_time_domain){
    if(! is.vector(time_domain))
      stop("'time_domain' variable must be a vector.")
    check.time_axis(time_domain)
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check coherence between flag for optimal parameters and optimal parameters
#'
#' @description
#' The function controls that one of the following condition is satisfied:
#' - if the flag for the optimal parameters is activated, then the optimal parameters should be provided in input
#' - if the flag is not activated, then the optimal parameters should not be provided and the parameter vector should be NU
#'
#' @param optimal_params Either a numerical vector of length equal to the number of model parameters or a NULL value.
#' @param flag_optimal_params Logical value. Did the user want to provide optimal parameters vector? If so, the variable should be TRUE; otherwise
#' (no optimal parameters), FALSE.
#'
#' @return An error if any condition is not satisfied.

check.flag_optimal_params <- function(optimal_params, flag_optimal_params){
  # If the flag is activated, then the user needs to provide the optimal parameter vector
  if(flag_optimal_params){
    if((!is.vector(optimal_params)) || (is.null(optimal_params)))
      stop("'optimal_params' must be provided with TRUE 'flag_optimal_params'.")
  }
  # If the flag is not activated, then no optimal parameters should be provided
  else{
    if((!is.null(optimal_params)) || (is.vector(optimal_params)))
      warning("'optimal_params' provided but FALSE 'flag_optimal_params'.")
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check structure for the Parameters Confidence Interval
#'
#' @description
#' The function controls that the structure of the Parameters Confidence Intervals coincides with the theoretical one.
#'
#' @param parametersCI S3 object of class 'ParametersCI', composed of two elements:
#' - left confidence interval: numerical vector of length equal to the number of parameters in the model
#' - right confidence interval: numerical vector of length equal to the number of parameters in the model
#'
#' @return An error if any condition is not satisfied.

check.structure_paramsCI <- function(parametersCI){
  # Control class
  if(!inherits(parametersCI, "ParametersCI"))
    stop("Wrong S3 class object for 'ParametersCI'.")
  
  # Control structure
  if(! is.list(parametersCI))
    stop("'ParametersCI' is not a list.")
  
  # Control its elements
  if(! is.vector(parametersCI$ParamsCI_left))
    stop("Wrong structure for 'left CI'.")
  
  if(! is.vector(parametersCI$ParamsCI_right))
    stop("Wrong structure for 'right CI'.")
}
#-------------------------------------------------------------------------------
check.params_range <- function(params_range, n_params){
  # Check structure of input class
  if(!inherits(params_range, "ParametersRange"))
    stop("'params_range' is not of class 'ParametersRange'.")
  
  # Check its number of elements
  if(length(params_range) != 2)
    stop("Wrong number of elements in object class 'ParametersRange'.")
  
  # Check dimensions and vector
  if(! is.vector(params_range$ParametersRangeMin))
    stop("'ParametersRangeMin' is not a vector.")
  if(! is.vector(params_range$ParametersRangeMax))
    stop("'ParametersRangeMiax' is not a vector.")
  
  if(length(params_range$ParametersRangeMin) != n_params)
    stop("Wrong length of 'ParametersRangeMin' vector.")
  if(length(params_range$ParametersRangeMax) != n_params)
    stop("Wrong length of 'ParametersRangeMax' vector.")
}
#-------------------------------------------------------------------------------
#' @title 
#' Check presence of null or nan element value in the dataset
#' 
#' @description
#' The method controls that the dataset does not contain 'NULL', 'null', 'NaN' or 'nan' value.
#' 
#' @param data Dataset (dataframe)
#'
#' @return An error if any condition is not satisfied.

check.dataset <- function(data){
  # Extract information from input variable
  n_col <- dim(data)[2]
  n_row <- dim(data)[1]
  
  # Check there are no null or nan value in the dataset
  for(j in 1:n_col){
    for(i in 1:n_row){
      if((data[i,j] == 'NaN') || (data[i,j] == 'nan') || 
         (data[i,j] == 'null') || (data[i,j] == 'NULL')){
        msg <- paste('Element in position (', i, ',',j,') is null or nan.')
        stop(msg)
      }
    }
  }
}
#-------------------------------------------------------------------------------
#' @title 
#' Check existence of provided input index
#' 
#' @description
#' The method controls that the provided input index exists: it cannot be greater than the maximum number
#' of parameters of the current model.
#' 
#' @param index Index with respect to which the user wants to study the one dimensional behaviour of 
#' the log-likelihood function.
#' @param n_params Number of parameters of the model
#'
#' @return An error if any condition is not satisfied.

check.index <- function(index, n_params){
  # Control existence of the input index
  if((index < 0) || (index > n_params)){
    stop('Provided index out of range.')
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Check positiveness of the multiplicative constant C
#' 
#' @description
#' The method controls the multiplicative constant C is non-negative (i.e. positive).
#' 
#' @param C_mult Multiplicative constant
#'
#' @return An error if the condition is not satisfied.

check.C_mult <- function(C_mult){
  # Controls the constant is positive
  if(C_mult <= 0)
    stop("Multiplicative constant is not-positive.")
}