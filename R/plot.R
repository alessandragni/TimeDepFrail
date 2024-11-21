
#-------------------------------------------------------------------------------
#' @title
#' Plot the Baseline Hazard Step-Function
#'
#' @description
#' This function plots the baseline hazard step-function based on the estimated parameters from the Adapted Paik et al.'s model.
#'
#' @details
#' The function plots a horizontal segment for each interval of the time domain, representing the baseline hazard. 
#' The boundaries of each segment are marked with colored dots, and subsequent segments are intentionally left unconnected 
#' to reflect the discrete nature of the intervals.
#'
#' @param result S3 object of class 'AdPaik', returned by the method call 'AdPaikModel(...)'.
#' @param xlim A numeric vector of length 2 specifying the x-axis limits. Default is set to 0 and the number of time-domain intervals.
#' @param ylim A numeric vector of length 2 specifying the y-axis limits. Default is 0 to the maximum value of the baseline hazard.
#' @param xlab,ylab String giving the x and y axis name. Default values are 'x' and 'y'.
#' @param main_title Title of the plot. Default title is 'Baseline hazard step-function'.
#' @param color Color used for plotting the horizontal segments of the step-function. Default one is 'black'.
#' @param pch Symbol for marking the boundaries of each segment. Default is a dot (value 21).
#' @param bg Color for the boundary symbols. Default matches the plot color ('black').
#' @param cex_points Size of the boundary symbols. Default is 0.7.
#'
#' @return Plot of the baseline hazard step-function and value of the function in each interval.
#'
#' @export
#' 
#' @examples
#' # Import data
#' data(data_dropout)
#' 
#' # Define the variables needed for the model execution
#' eps_paik <- 1e-10
#' categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
#' categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#'
#'\donttest{
#' # Call the main model function
#' result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#'
#' plot_bas_hazard(result)
#' }
plot_bas_hazard <- function(result,
                           xlim = c(0,length(result$TimeDomain)-1), ylim = c(0,max(result$BaselineHazard)),
                           xlab = "x", ylab = "y", main_title = "Baseline hazard step-function",
                           color = "black", pch = 21, bg = "black", cex_points = 0.7){

  # Check correctness of result structure
  check.result(result)

  # Extract information from input variables
  time_axis <- result$TimeDomain
  L <- n_intervals <- result$NIntervals
  optimal_params <- result$OptimalParameters
  eps <- 1e-2

  # Compute the baseline hazard function
  baseline_hazard <- result$BaselineHazard
  
  # Check it is non negative
  for(k in 1:L){
    if(baseline_hazard[k] < 0){
      msg <- paste("Negative baseline hazard value in position ", k, ".")
      stop(msg)
    }
  }

  # Plot the baseline hazard using a horizontal segment for each interval
  # dev.new()
  plot(c(time_axis[1], time_axis[2] - eps), c(baseline_hazard[1], baseline_hazard[1]), col = color,
       main = main_title, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       cex = cex_points, pch = pch, bg = "black")
  lines(c(time_axis[1], time_axis[2]-eps), c(baseline_hazard[1], baseline_hazard[1]), col= color)
  for(i in 2:L){
    points(time_axis[i] + eps, baseline_hazard[i], col = color, cex = cex_points, pch = pch, bg = "black")
    points(time_axis[i+1] - eps, baseline_hazard[i], col = color, cex = cex_points, pch = pch, bg = "black")
    lines(c(time_axis[i] + eps, time_axis[i+1] - eps), c(baseline_hazard[i], baseline_hazard[i]), col = color)
  }
}
#-------------------------------------------------------------------------------
#' @title
#' Plot the Posterior Frailty Estimates
#'
#' @description
#' This function plots the posterior frailty estimates for each group in each time interval. 
#' Each group’s estimates are represented by a sequence of points connected by straight lines. 
#' The function can plot either the entire posterior frailty estimate or 
#' its time-independent and time-dependent components based on user-specified flags.
#'
#' @details
#' Recalling the frailty structure as \eqn{Z_{jk} = \alpha_{j} + \epsilon_{jk}, \forall j,k} and the posterior
#' frailty estimate as \eqn{\hat{Z}_{jk} = \hat{\alpha}_{j}/\hat{\alpha}_{max} + \hat{\epsilon}_{jk}/\hat{\epsilon}_{max}}, 
#' this function allows plotting either the entire posterior frailty estimate \eqn{\hat{Z}_{jk}}
#' or its time-independent \eqn{\frac{\hat{\alpha}_{j}}{\hat{\alpha}_{\text{max}}}} or 
#' time-dependent \eqn{\frac{\hat{\epsilon}_{jk}}{\hat{\epsilon}_{\text{max}}}} components.
#' The user can control which components to display using the flag_eps and flag_alpha parameters. 
#' Only one of these flags can be set to TRUE at a time.
#'
#' @param result S3 object of class 'AdPaik', returned by the method call 'AdPaikModel(...)'.
#' @param data Dataset (dataframe) in which all variables of the formula object must be found and contained.
#' @param flag_eps Logical flag indicating whether to plot only the time-dependent posterior frailty estimates. Default is FALSE.
#' @param flag_alpha Logical flag indicating whether to plot only the time-independent posterior frailty estimates. Default is FALSE.
#' @param xlim Numeric vector of length 2 specifying the x-axis limits. Default is based on the number of intervals.
#' @param ylim Numeric vector of length 2 specifying the y-axis limits. Default is (0, 10).
#' @param xlab,ylab String giving the x and y axis name. Default values are 'Intervals' and 'Values'.
#' @param main_title Title of the plot. Default title is 'Posterior frailty estimates'.
#' @param cex Dimension of the points used for plotting the estimates.
#' @param pch_type Numerical vector of length equal to the number of clusters in the data, giving the symbol to be used for plotting the estimates.
#' Default symbol (circle, 21) is the same for all clusters.
#' @param color_bg Numerical vector of length equal to the number of clusters in the data, giving the color to be used for plotting the symbols
#' for the estimates. Default ('black') is the same for all faculties. On the other hand, the same color is used throughout the intervals for
#' the same faculty.
#' @param cex_legend Dimension of the symbol in the legend. Default is 0.7.
#' @param pos_legend  Either a numeric vector providing the x and y coordinates for the legend or 
#' a string specifying the legend's position (e.g., 'bottomright', 'bottom', 'bottomleft', 'left',
#' 'topleft', 'top', 'topright', 'right', 'center').
#'
#' @return The plot of the posterior frailty estimates.
#'
#' @export
#' 
#' @examples
#' # Import data
#' data(data_dropout)
#' 
#' # Define the variables needed for the model execution
#' eps_paik <- 1e-10
#' categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
#' categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#'
#' # Call the main model function
#' 
#' \donttest{
#' result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#'
#' # Define variables for plotting the estimates
#' pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
#' color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
#' 
#' plot_post_frailty_est(result, data_dropout,
#'                       xlim = c(1, 11), ylim = c(0,3), 
#'                       pch_type = pch_type, color_bg = color_bg,
#'                       xlab = 'Time [intervals]', ylab = 'Posterior estimates',
#'                       pos_legend = 'bottomright')
#'  }                     

plot_post_frailty_est <- function(result, data,
                                  flag_eps = FALSE, flag_alpha = FALSE,
                                  xlim = c(0,length(result$TimeDomain)-1), ylim = c(0, 10),
                                  xlab = "Intervals", ylab = "Values", main_title = "Posterior frailty estimates",
                                  cex = 0.7,
                                  pch_type = rep(21, length(centre_codes)),
                                  color_bg = rep("black", length(centre_codes)),
                                  cex_legend = 0.7, pos_legend = "topright"){

  # Check correctness of result structure
  check.result(result)

  # Extract information from input variables
  time_axis <- result$TimeDomain
  L <- n_intervals <- result$NIntervals
  formula <- result$formula
  post_frailty_est <- result$PosteriorFrailtyEstimates

  # Check variables of formula are contained in dataset
  check.formula_terms(formula, data)

  # Extract elements of the formula object
  formula_variables <- all.vars(formula)

  # Extract position of cluster
  special <- c("cluster")
  terms_object <- terms(formula, specials = special, data = data)
  cluster_index <- attr(terms_object, "specials")$cluster
  cluster_name <- formula_variables[cluster_index]

  # Extract centre_codes
  centre <- data[,cluster_name]
  centre_codes <- levels(factor(centre))
  n_centres <- length(centre_codes)

  # Check correctness of post_frailty_list and centre
  check.structure_post_frailty_est(post_frailty_est, n_intervals, n_centres)
  check.value_post_frailty(post_frailty_est, n_centres, n_intervals)
  check.centre(centre_codes)

  # Control that at most one between flag_eps and flag_alpha is TRUE
  if((flag_eps == TRUE) & (flag_alpha == TRUE))
    stop("At most one flag must be TRUE, either 'flag_eps' or 'flag_alpha'")

  # Check correctness of pch_type and color_bg variables
  check.pchtype_colorbg(centre_codes, pch_type, color_bg)

  # Check correctness of pos_legend
  check.poslegend(pos_legend)

  # Define what to plot, according to the flag
  post_fralty <- 0
  if(flag_eps == TRUE)
    post_frailty <- post_frailty_est$eps
  if(flag_alpha == TRUE)
    post_frailty <- post_frailty_est$alpha
  if((flag_eps == FALSE) & (flag_alpha == FALSE))
    post_frailty <- post_frailty_est$Z

  # Plot posterior frailty estimates
  # dev.new()
  if(flag_alpha == FALSE) {
    plot(seq(1,n_intervals,1), post_frailty[1,],
         pch = pch_type[1], bg = color_bg[1], cex = cex,
         main = main_title, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim)
    for(k in 1:(n_intervals-1))
      lines(c(k,k+1), c(post_frailty[1,k], post_frailty[1,k+1]))
    for(i in 2:n_centres){
      for(k in 1:(n_intervals-1)){
        points(k, post_frailty[i,k], pch = pch_type[i], bg = color_bg[i], cex = cex)
        points(k+1, post_frailty[i,k+1], pch = pch_type[i], bg = color_bg[i], cex = cex)
        lines(c(k,k+1), c(post_frailty[i,k], post_frailty[i,k+1]))
      }
    }
  } else {
    plot(rep(1, length(post_frailty)), post_frailty,
         pch = pch_type[1], bg = color_bg[1], cex = cex,
         main = main_title, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim)
    for(i in 2:n_centres){
        points(1, post_frailty[i], pch = pch_type[i], bg = color_bg[i], cex = cex)
    }
  }
  
  if(is.vector(pos_legend))
    legend(pos_legend[1], pos_legend[2], legend = centre_codes, col = color_bg,
           pch = pch_type, pt.bg = color_bg, cex = cex_legend)
  if(is.character(pos_legend))
    legend(pos_legend, legend = centre_codes, col = color_bg,
           pch = pch_type, pt.bg = color_bg, cex = cex_legend)
}

#-------------------------------------------------------------------------------
#' @title
#' Plot for the Frailty Standard Deviation or Variance
#'
#' @description
#' This function generates a plot of either the frailty standard deviation or the frailty variance for the intervals in the time-domain.
#'
#' @details
#' The plot represents the values of the frailty standard deviation or variance for each time interval. 
#' It connects these points to illustrate the trend of the chosen metric.
#'
#' This function supports two modes of operation:
#' - Plotting the frailty standard deviation or variance retrieved from the main model (contained in the S3 object of class 'AdPaik').
#' - Plotting a user-provided vector of frailty standard deviations, which can be computed using the method `frailty.sd`. This allows for flexibility in analysis without re-optimizing the log-likelihood function. For instance, users can compare frailty standard deviations computed with different model specifications (e.g., including only time-dependent terms).
#'
#' The output will differentiate between these two cases, ensuring the correct values are plotted regardless of the source.
#'
#'
#' @param result An S3 object of class 'AdPaik', returned by the main model call 'AdPaikModel(...)'.
#' @param frailty_sd A numerical vector representing the evaluated frailty standard deviation, with length equal to the number of time-domain intervals.
#' Its elements must be non-negative.
#' @param flag_variance A boolean flag indicating whether to plot the frailty variance (`TRUE`) or the frailty standard deviation (`FALSE`). Default is `FALSE`.
#' @param flag_sd_external A logical flag indicating whether the user is providing an external frailty standard deviation vector.
#' @param xlim A numeric vector of length 2, defining the range for the x-axis (intervals). Default is from 1 to the number of intervals.
#' @param ylim A numeric vector of length 2, defining the range for the y-axis (values). Default is `(0, 10)`.
#' @param xlab A string for the x-axis label. Default is `'Intervals'`.
#' @param ylab A string for the y-axis label. Default is `'Values'`.
#' @param main_title A string for the plot title. Default title is `'Frailty Standard Deviation'`.
#' @param pch A numeric or character symbol used for plotting the frailty standard deviation values. Default is a dot (`21`).
#' @param color_bg A string specifying the color used for the symbols. Default is `'blue'`.
#' @param cex_points A numeric value indicating the size of the plotting symbols. Default is `0.7`.
#'
#' @return A plot displaying either the frailty standard deviation or variance across the specified intervals.
#' 
#' @export
#' 
#' @examples
#' # Import data
#' data(data_dropout)
#' 
#' # Define the variables needed for the model execution
#' eps_paik <- 1e-10
#' categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
#' categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#'
#' # Call the main model function
#' \donttest{
#' result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#'
#' plot_frailty_sd(result, ylim=c(0, 0.50), xlab = 'Time [intervals]', ylab = 'Standard deviation')
#' }

plot_frailty_sd <- function(result, frailty_sd = NULL, flag_variance = FALSE, flag_sd_external = FALSE,
                            xlim = c(1, length(result$TimeDomain)-1), ylim = c(0, 10),
                            xlab = "Intervals", ylab = "Values", main_title = "Frailty standard deviation",
                            pch = 21, color_bg = "blue", cex_points = 0.7){
  # Check result
  check.result(result)

  # Extract information from input variables and initialize frailty standard deviation
  L <- n_intervals <- result$NIntervals
  values <- c()

  # Check coherence between input variables
  if(flag_sd_external){
    if(is.null(frailty_sd))
      stop("Expected standard deviation vector provided by the user.")
    else{
      check.pos_frailty_sd(frailty_sd, n_intervals)
    }
  }

  # Initialize values vector
  if(! flag_sd_external){
    if(flag_variance)
      values <- result$FrailtyDispersion$FrailtyVariance
    else
      values <- result$FrailtyDispersion$FrailtyStandardDeviation
  }
  else{
    if(flag_variance)
      values <- frailty_sd
    else
      values <- (frailty_sd)^2
  }

  # Plot standard deviation of the frailty
  # dev.new()
  plot(1, values[1], pch = pch, bg = color_bg,
       xlab = xlab, ylab = ylab, main = main_title,
       xlim = xlim, ylim = ylim)
  for (k in 2:L){
    points(k, values[k], pch = pch, bg = color_bg, cex = cex_points)
    lines(c(k-1,k), c(values[k-1], values[k]))
  }
}

#-------------------------------------------------------------------------------
#' @title
#' Plot the One-Dimensional Log-Likelihood Function
#'
#' @description
#' This function plots the trend of the log-likelihood function with respect to a single parameter specified by its index in the parameter vector. 
#' It generates samples of the parameter, evaluates them in the log-likelihood function, and displays the results along with the maximum point of the one-dimensional log-likelihood function.
#'
#' @param param_1D A numeric value representing the optimal parameter determined by maximizing the log-likelihood function for the specified parameter.
#' @param index_param_1D An integer representing the index of the optimal parameter within the parameter vector.
#' @param ll_1D A numeric value of the log-likelihood function evaluated at the optimal parameter `param_1D`, with other parameters held constant.
#' @param params A numeric vector of length equal to the number of parameters minus one, containing fixed values for the other parameters.
#' @param param_range_min A numeric value indicating the minimum allowable value for the parameter `param_1D`.
#' @param param_range_max A numeric value indicating the maximum allowable value for the parameter `param_1D`.
#' @param dataset A data frame or matrix containing individual covariates.
#' @param centre A numeric vector indicating individual cluster membership; its length must match the number of individuals in the dataset.
#' @param time_axis A numeric vector corresponding to the subdivisions of the temporal domain.
#' @param dropout_matrix A binary matrix indicating which interval of the time domain an individual failed. Each row should sum to 1 (if failed) or 0 (if not failed), with dimensions (n_individuals, n_intervals).
#' @param e_matrix A matrix of dimensions (n_individuals, n_intervals) where each element contains the evaluation of the temporal integral performed by the function `time_int_eval`.
#' @param n_points An integer specifying the number of points at which to evaluate the log-likelihood function. A value that is neither too small nor too high is recommended; the default is 150.
#' @param cex A numeric value specifying the size of the points used for the graphical representation of the log-likelihood function. Default is 0.7.
#' @param cex_max A numeric value indicating the size of the optimal point (the one maximizing the log-likelihood function). Default is 0.8.
#' @param color_bg A string specifying the color for the points representing the log-likelihood trend. Default is `'black'`.
#' @param color_max_bg A string specifying the color for the optimal point provided as the first argument. Default is `'red'`.
#' @param pch A numeric or character symbol representing the shape of the plotted points. Default is a circle (`21`).
#'
#' @return A plot displaying the trend of the log-likelihood function with respect to a single parameter, including the maximum point.
#'
#' @export
plot_ll_1D.AdPaik <- function(param_1D, index_param_1D, ll_1D, params, param_range_min, param_range_max,
                              dataset, centre, time_axis, dropout_matrix, e_matrix,
                              n_points = 150,
                              cex = 0.7, cex_max = 0.8, color_bg = "black", color_max_bg = "red",
                              pch = 21){
  
  # Define the structure containing the generated points and the associated log-likelihood value
  #param_values <- rep(0, n_points)
  ll_values <- rep(0, n_points)
  
  # Generate n_points for the indicated parameter inside its min, max range
  param_values <- runif(n_points, param_range_min, param_range_max)
  
  # For each point, evaluate the log-likelihood function
  for(i in 1:n_points){
    params[index_param_1D] <- param_values[i]
    ll_values[i] <- ll_AdPaik_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix)
  }
  
  # Plot the log-likelihood trend with respect to the indicated parameter
  string_title <- paste("Log-likelihood trend wrt parameter ", index_param_1D)
  
  # dev.new()
  plot(param_values, ll_values, pch=pch, col=color_bg, cex = cex,
       xlim = c(param_range_min, param_range_max), #ylim=c(min(ll_values), max(ll_values)),
       main = string_title, xlab = "Values", ylab = "Log-likelihood")
  points(param_1D, ll_1D, bg = color_max_bg, pch = pch, cex = cex_max)
}



#-------------------------------------------------------------------------------
#' @title
#' Plot the One-Dimensional Log-Likelihood Function
#'
#' @description
#' This function plots the trend of the log-likelihood function concerning a single parameter specified by its index in the parameter vector. 
#' It generates samples of the parameter, evaluates them in the log-likelihood function, and displays the results along with the maximum point of the one-dimensional log-likelihood function.
#'
#' @param param_1D A numeric value representing the optimal parameter determined by maximizing the log-likelihood function for the specified parameter.
#' @param index_param_1D An integer representing the index of the optimal parameter within the parameter vector.
#' @param ll_1D A numeric value of the log-likelihood function evaluated at the optimal parameter `param_1D`, with the other parameters held constant.
#' @param params A numeric vector of length equal to the number of parameters minus one, containing the fixed values for the other parameters.
#' @param param_range_min A numeric value indicating the minimum allowable value for the parameter `param_1D`.
#' @param param_range_max A numeric value indicating the maximum allowable value for the parameter `param_1D`.
#' @param dataset A data frame or matrix containing individual covariates.
#' @param centre A numeric vector indicating individual cluster membership; its length must match the number of individuals in the dataset.
#' @param time_axis A numeric vector corresponding to the subdivisions of the temporal domain.
#' @param dropout_matrix A binary matrix indicating which interval of the time domain an individual failed. Each row should sum to 1 (if failed) or 0 (if not failed), with dimensions (n_individuals, n_intervals).
#' @param e_matrix A matrix of dimensions (n_individuals, n_intervals), where each element contains the evaluation of the temporal integral performed by the function `time_int_eval`.
#' @param n_points An integer specifying the number of points at which to evaluate the log-likelihood function. A value that is neither too small nor too high is recommended; the default is 150.
#' @param cex A numeric value specifying the size of the points used for the graphical representation of the log-likelihood function. Default is 0.7.
#' @param cex_max A numeric value indicating the size of the optimal point (the one maximizing the log-likelihood function). Default is 0.8.
#' @param color_bg A string specifying the color for the points representing the log-likelihood trend. Default is `'black'`.
#' @param color_max_bg A string specifying the color for the optimal point provided as the first argument. Default is `'red'`.
#' @param pch A numeric or character symbol representing the shape of the plotted points. Default is a circle (`21`).
#'
#' @return A plot displaying the trend of the log-likelihood function concerning a single parameter, including the maximum point.
#'
#' @export
plot_ll_1D <- function(param_1D, index_param_1D, ll_1D, params, param_range_min, param_range_max,
                       dataset, centre, time_axis, dropout_matrix, e_matrix,
                       n_points = 150,
                       cex = 0.7, cex_max = 0.8, color_bg = "black", color_max_bg = "red",
                       pch = 21){
  
  # Define the structure containing the generated points and the associated log-likelihood value
  #param_values <- rep(0, n_points)
  ll_values <- rep(0, n_points)
  
  # Generate n_points for the indicated parameter inside its min, max range
  param_values <- runif(n_points, param_range_min, param_range_max)
  
  # For each point, evaluate the log-likelihood function
  for(i in 1:n_points){
    params[index_param_1D] <- param_values[i]
    ll_values[i] <- ll_AdPaik_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix)
  }
  
  # Plot the log-likelihood trend with respect to the indicated parameter
  string_title <- paste("Log-likelihood trend wrt parameter ", index_param_1D)
  
  # dev.new()
  plot(param_values, ll_values, pch=pch, col=color_bg, cex = cex,
       xlim = c(param_range_min, param_range_max), #ylim=c(min(ll_values), max(ll_values)),
       main = string_title, xlab = "Values", ylab = "Log-likelihood")
  points(param_1D, ll_1D, bg = color_max_bg, pch = pch, cex = cex_max)
}
#-------------------------------------------------------------------------------

