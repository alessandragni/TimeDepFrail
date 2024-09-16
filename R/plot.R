#' Function for plotting the baseline hazard step-function, using the estimated parameters.
#' For each interval of the time-domain a horizontal segment is traced and subsequent elements
#' are connected through a line (default dotted), which dimension and type are user-specified.
#' To have a nice graphical representation, and not a vertical line connecting two subsequent steps,
#' it is suggested to shrink the internal segments of a small quantity (@eps), on both sides.
#'
#' @param result S3 object of class 'AdPaik', returned by the method call 'AdPaikModel(...)'.
#' @param xlim Limits of the horizontal axis (x)
#' @param ylim Limits of the vertical axis (y)
#' @param xlab Name for the horizontal axis (x)
#' @param ylba Name for the vertical axis (y)
#' @param main_title Title for the plot
#' @param color Color to be used for the graphical representation of the baseline hazard step-function
#' @param name descriptioncex_points Dimension of the points, delimiting each interval horizontal segment of the plot
#' @param name descriptionlty Type of line connecting two consecutive segments (default (3) dotted)
#' @param name descriptionlwd Width of the line connecting two consecutive segments (default 1)
#' @param eps Small quantity to be subtracted on both sides of the internal segments. Default value equal to 5e-2.

# plot.bas_hazard_con <- function(result,
#                             xlim = c(0,time_axis[length(time_axis)]), ylim = c(0,1),
#                             xlab = "x", ylab = "y", main_title = "Baseline hazard step-function",
#                             color = "black", pch = 21, bg = "black", cex_points = 0.7,
#                             lty = 3, lwd = 1, eps = 5e-2){
#
#   # Check correctness of result structure
#   check.result(result)
#
#   # Extract information from input variables
#   time_axis <- result$TimeDomain
#   L <- n_intervals <- result$NIntervals
#   optimal_params <- result$OptimalParameters
#
#   # Compute the baseline hazard function
#   baseline_hazard <- result$BaselineHazard
#
#   # Plot the baseline hazard using a horizontal segment for each interval
#   dev.new()
#   plot(c(time_axis[1], time_axis[2] - eps), c(baseline_hazard[1], baseline_hazard[1]), col = color,
#        main = main_title, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
#        cex = cex_points, pch = pch, bg = "black")
#   lines(c(time_axis[1], time_axis[2]-eps), c(baseline_hazard[1], baseline_hazard[1]), col= color)
#   for(i in 2:L){
#     points(time_axis[i] + eps, baseline_hazard[i], col = color, cex = cex_points, pch = pch, bg = "black")
#     points(time_axis[i+1] - eps, baseline_hazard[i], col = color, cex = cex_points, pch = pch, bg = "black")
#     lines(c(time_axis[i] + eps, time_axis[i+1] - eps), c(baseline_hazard[i], baseline_hazard[i]), col = color)
#   }
#   for(i in 2:L){
#     lines(c(time_axis[i]-eps, time_axis[i]+eps), c(baseline_hazard[i-1], baseline_hazard[i]),
#           lty = lty, lwd = lwd, col = color)
#   }
#   points(time_axis[1], 0.0, col = color, cex = cex_points, pch = pch, bg = "black")
#   lines(c(time_axis[1], time_axis[1]), c(0.0, baseline_hazard[1]),
#         lty = lty, lwd = lwd, col = color)
# }

#-------------------------------------------------------------------------------
#' @title
#' Plot the baseline hazard step-function
#'
#' @description
#' Function for plotting the baseline hazard step-function, using the estimated parameters.
#'
#' @details
#' For each interval of the time-domain a horizontal segment is traced and its boundaries are marked with a full colored dot.
#' Subsequent segments are not connected, for a design matter.
#'
#' @param result S3 object of class 'AdPaik', returned by the method call 'AdPaikModel(...)'.
#' @param xlim Numeric vector of length 2, giving the x coordinate ranges.
#' Default value goes from 0 up to the number of intervals of the time-domain.
#' @param ylim Numeric vector of length 2, giving the y coordinate ranges.
#' Deafult value goes from 0 up to the maximum value reached by the baseline hazard.
#' @param xlab,ylab String giving the x and y axis name.
#' Deafult values are 'x' and 'y'.
#' @param main_title Title of the plot. Default title is 'Baseline hazard step-function'.
#' @param color Color used for plotting the horizontal segments of the step-function. Deafult one is 'black'.
#' @param pch Symbol for marking the bounaries of each segment. Deafult one is a dot (number 21).
#' @param bg Color used for plotting the boundary symbol. Default one is equal to the plot color ('black').
#' @param cex_points Dimension of the boundary symbol. Deafult is 0.7.
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
#' # Call the main model function
#' result <- AdPaikModel(formula, data, time_axis, categories_range_min, categories_range_max)
#'
#' plot_bas_hazard(result)
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
#' Plot the posterior frailty estimates
#'
#' @description
#' Function for plotting the posterior frailty estimate of each group in each interval,
#' using a single point. For the same group, the sequence of points is then connected using a straight
#' line.
#'
#' @details
#' Recalling the frailty structure as \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k} and the posterior
#' frailty estimate as \eqn{\hat{Z}_{jk} = \hat{\alpha}_j/\hat{\alpha}_{max} + \hat{\epsilon}_{jk}/\hat{\epsilon}_{max}}, it is
#' possible to plot the entire posterior frailty estimate \eqn{\hat{Z}_{jk}} or just a part: either the time-dependent estimate
#' \hat{\epsilon}_{jk}/\hat{\epsilon}_{max} or the time-independent one \hat{\alpha}_j/\hat{\alpha}_{max}.
#' To do so, the user has to specify the alue of two flags.

#' @param result S3 object of class 'AdPaik', returned by the method call 'AdPaikModel(...)'.
#' @param data Dataset (dataframe) in which all variables of the formula object must be found and contained.
#' @param flag_eps Does the user want to plot only the time-dependent posterior frailty estimates? If so, the flag must be TRUE.
#' Otherwise, for the entire posterior estimates, it must be FALSE (default value).
#' @param flag_alpha Does the user want to plot only the time-independent posteior frailty estimates? If so, the flag must be TRUE.
#' Otherwise, for the entire posterior estimates, it must be FALSE (default value).
#' Note that, it is not possible to have both previous flags TRUE: either one of the two must be TRUE. However, both can be FALSE.
#' @param xlim Numeric vector of length 2, giving the x coordinate ranges.
#' Default value goes from 1 up to the number of intervals of the time-domain, because a point is plotted for each faculty in each interval.
#' @param ylim Numeric vector of length 2, giving the y coordinate ranges.
#' Default value is (0,10). Pay attention that no non-negative values only can be accepted.
#' @param xlab,ylab String giving the x and y axis name. Default values are 'Intervals' and 'Values'.
#' @param main_title Title of the plot. Default title is 'Posterior frailty estimates'.
#' @param cex Dimension of the points used for plotting the estimates.
#' @param pch_type Numerical vector of length equal to the number of clusters in the data, giving the symbol to be used for plotting the estimates.
#' Default symbol (circle, 21) is the same for all clusters.
#' @param color_bg Numerical vector of length equal to the number of clusters in the data, giving the color to be used for plotting the symbols
#' for the estimates. Default ('black') is the same for all faculties. On the other hand, the same color is used throughout the intervals for
#' the same faculty.
#' @param cex_legend Dimension of the symbol in the legend. Default is 0.7.
#' @param pos_legend It could be either a numerical vector of length 2, providing the x and y coordinates of the legend, or a string indicating where
#' the legend has to be located in the plot. In the last case, the only possibilities are: 'bottomright', 'bottom', 'bottomleft', 'left',
#' 'topleft', 'top', 'topright', 'right', 'center'.
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
#' result <- AdPaikModel(formula, data, time_axis, categories_range_min, categories_range_max)
#'
#' # Define variables for plotting the estimates
#' pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
#' color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
#' 
#' plot_post_frailty_est(result, data_dropout, ylim=c(0,2),
#'                       pch_type = pch_type, color_bg = color_bg)
plot_post_frailty_est <- function(result, data,
                                  flag_eps = FALSE, flag_alpha = FALSE,
                                  xlim = c(1, length(time_axis) - 1), ylim = c(0, 10),
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
  terms_object <- terms(formula, special, data = data)
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
#' Plot for the frailty standard deviation
#'
#' @description
#' Function for plotting either the frailty standard deviation or the frailty variance, according to the value
#' of a precise and user-provided flag.
#'
#' @details
#' For each interval of the time-domain, a single point is plotted with value at the computed frailty standard deviation or variance.
#' Subsequent points are connected to show the trend of the required function.
#'
#' @details
#' The method gives the possibility of plotting either
#' - the frailty standard deviation returned by the main model call and, therefore, contained in the S3 object class 'AdPaik'
#' or
#' - the frailty standard deviation computed through the method 'frailty.sd'. We recall this method is introduced to permit the user to change the
#' way the frailty standard deviation has been computed, without the necessity of optimizing again the log-likelihood function.
#' For instance, we require the computation of the frailty standard deviation with both terms (time-dependent and time-independent) in the main model
#' 'AdPaikModel', but now we want also the frailty standard deviation with the sole time-dependent term. Thanks to the introduced method 'frailty.sd',
#' we can compute again the frailty standard deviation without maximizing again the log-likelihood function.
#' The result we obtain in the second case is different in values and structure, because we simply get a vector of L elements instead
#' of a class, and we want the plot to be able to distinguish the case and plot it anyway.
#'
#'
#' @param result S3 object of class 'AdPaik', returned by the main model call 'AdPaikModel(...)'.
#' @param frailty_sd Numerical vector of evaluated frailty standard deviation, of length equal to the number of intervals of the time-domain.
#' Its element must be non-negative.
#' @param flag_variance Do you want to plot the frailty standard deviation? If so, the flag must be equal to TRUE; otherwise, to FALSE (if we want to
#' plot the variance). Default value is FALSE.
#' @param flag_sd_external Logical value. Do you want to provide a frailty standard deviation
#' @param xlim Numerical vector of length 2, giving the x coordinate ranges.
#' Default value goes from 1 up to the number of intervals of the time-domain, because for each interval a single point is plotted and not a segment.
#' @param ylim Numerical vector of length 2, giving the y coordinate ranges.
#' Default value is (0,10).
#' @param xlab,ylab String giving the x and y axis name.
#' Deafult values are 'Intervals' and 'Values'.
#' @param main_title Title of the plot. Default title is 'Frailty standard deviation'.
#' @param pch Symbol used for plotting the frailty standard deviation value. Default is a dot (21).
#' The same symbol is used throughout the plot, with no possibility of change.
#' @param color_bg Color used for the symbols. Deafult is 'blue'.
#' @param cex_points Dimension of the symbols. Deafult is 0.7.
#'
#' @return Plot of either the frailty standard deviation or the frailty variance.
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
#' result <- AdPaikModel(formula, data, time_axis, categories_range_min, categories_range_max)
#'
#' plot_frailty_sd(sd, time_axis, FALSE, ylim = c(0,1))
plot_frailty_sd <- function(result, frailty_sd = NULL, flag_variance = FALSE, flag_sd_external = FALSE,
                            xlim = c(1, length(time_axis)-1), ylim = c(0, 10),
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
}#-------------------------------------------------------------------------------
#' @title
#' Plot the one-dimensional log-likelihood function
#'
#' @description
#' Function for plotting the trend of the log-likelihood function with respect to
#' a single parameter, whose position in the parameter vector is specified in the argument call.
#' For such plot, some samples of the parameter under consideration are generated and
#' then evaluated in the log-likelihood function. These are plotted together with the maximum point of the
#' one-dimensional log-likelihood function, provided in input.
#'
#' @param param_1D Optimal parameter value determined maximizing the log-likelihood function with respect to it.
#' @param index_param_1D Index of the optimal parameter inside the parameter vector.
#' @param ll_1D Log-likelihood function evaluated at the optimal parameter @param_1D, with the other parameters assuming a fixed value.
#' @param params Numerical vector of length equal to the number of parameters minus one. It contains the fixed value of the other parameters.
#' @param param_range_min Minimum value assumable by the parameter param_1D.
#' @param param_range_max Maximum value assumable by the parameter param_1D.
#' @param dataset Dataset with individual covariates. It can be either a matrix or a dataframe.
#' @param centre Numerical vector of individual cluster membership. It must be length as the number of individual in the dataset.
#' @param time_axis Numerical vector of length 2, corresponding to the subdivision of the temporal domain.
#' @param dropout_matrix Binary matrix indicating in which interval of the time domain an individual failed. For each individual,
#' the sum of the row elements must be equal to 1 (if he/she failed) or 0 (if he/she does not failed).
#' It has dimension equal to (n_individuals, n_intervals)
#' @param e_matrix Matrix of dimension (n_individual, n_intervals), where each element contains the evaluation of the temporal
#' integral, performed through the function @time_int_eval.
#' @param n_points Number of points in which the log-likelihood function must be evaluated and then plotted.
#' To have a nice graphical representation, chose an intermediate value: not to small and not too high. Default value is 150.
#' @param cex Dimension of the points used for the graphical representation of the log-likelihood function. Default value is 0.7.
#' @param cex_max Dimension of the optimal point provided as first argument (i.e. the point that maximizes the log-likelihood
#' function from the point of view of a single parameter). Default value is 0.8.
#' @param color_bg Color to be used for the points of the log-likelihood trend. Deafult is 'black'.
#' @param color_max_bg Color to be used for the optimal point provided as first argument. Default is 'red'.
#' @param pch Shape of the plotted point. Deafult is a circle (21).
#'
#' @return Plot of the log-likelihood trend from the point of view of a single parameter.
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
       xlim = c(param_range_min, param_range_max), #, ylim=c(min(ll_values), max(ll_values)),
       main = string_title, xlab = "Values", ylab = "Log-likelihood")
  points(param_1D, ll_1D, bg = color_max_bg, pch = pch, cex = cex_max)
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
