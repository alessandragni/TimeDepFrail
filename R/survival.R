#' @title
#' Compute the Conditional Survival Function
#'
#' @description
#' Computes the conditional survival function based on the 'Adapted Paik et al.' model's 
#' given the estimated coefficients and frailty effects.
#'
#' @param result S3 object of class 'AdPaik' containing model results.
#'
#' @return A dataset where each row corresponds to an individual unit in the dataset, 
#' and the columns represent the survival function values over time interval, 
#' with the first column indicating the cluster to which the individual belongs. 
#'
#' @export
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
#' survivalAdPaik(result)
#'  }  
survivalAdPaik <- function(result) {
  # Check for valid input
  if (!inherits(result, "AdPaik")) stop("'result' must be of class 'AdPaik'.")
  
  data = data.frame(result$dataset)
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  
  # Extract beta coefficients and formula
  beta <- coef.AdPaik(result)  # Extract coefficients
  # names(coef(result)$beta) # Extract names
  full_formula <- result$formula  # Extract the full formula
  
  # Extract covariates excluding cluster terms
  terms_object <- terms(full_formula)
  covariates <- attr(terms_object, "term.labels")
  covariates_nocluster <- covariates[!grepl("^cluster\\(", covariates)]
  
  # Construct a design matrix for covariates, accounting for factors/dummies
  covariate_data <- data[, covariates_nocluster, drop = FALSE]
  design_matrix <- stats::model.matrix(~ . - 1, data = covariate_data)  # No intercept
  
  # Match design matrix columns with beta names
  if (is.null(names(beta))) {
    stop("Beta coefficients must have names to align with design matrix columns.")
  }
  
  # Select only the design matrix columns that match the beta names
  matching_columns <- colnames(design_matrix) %in% names(beta)
  if (!any(matching_columns)) {
    stop("No matching columns found between design matrix and beta coefficients.")
  }
  design_matrix <- design_matrix[, matching_columns, drop = FALSE]

  # Compute the linear predictor and exponentiate
  linear_predictor <- as.matrix(design_matrix) %*% beta
  exp_linear_predictor <- exp(c(linear_predictor))
  
  # Compute the the inner terms of the integral
  exp_phi <- exp(result$OptimalParameters[1:result$NIntervals])  # Baseline hazard parameters
  time_diffs <- diff(result$TimeDomain)  # Differences in the time domain
  posterior_frailty <- result$PosteriorFrailtyEstimates$Z  # Posterior frailty estimates
  inner_part <- t(t(posterior_frailty) * (exp_phi * time_diffs))  # Scale by hazard params
  
  # Compute cumulative hazard
  comput <- t(apply( inner_part , 1, cumsum)) 
  
  df_explinerarpred = data.frame('group'=data[[result$ClusterVariable]], 
                              'exp_linear_predictor'=exp_linear_predictor)
  df_comput = data.frame('group'=levels(factor(data[[result$ClusterVariable]])), 
                         'survival'=comput)
  
  # Perform the merge
  df_explinerarpred$row_id <- seq_len(nrow(df_explinerarpred))  # Add a row identifier
  result_df <- merge(df_explinerarpred, df_comput, by = "group", all.x = TRUE, sort = FALSE)
  # Restore the original order
  result_df <- result_df[order(result_df$row_id), ]
  result_df$row_id <- NULL  # Remove the temporary identifier if not needed
  
  survival = data.frame('group'=result_df$group,
                        t(apply(result_df$exp_linear_predictor * result_df[,3:ncol(result_df)], 1, function(row) exp(-row))))
  rownames(survival) = seq_len(nrow(survival))
  # Return survival function
  return(survival)
    #as.matrix(survival[,2:ncol(survival)]))
}


#-------------------------------------------------------------------------------
#' @title
#' Plot of Conditional Survival Function
#'
#' @description
#' Plots the conditional survival function based on the 'Adapted Paik et al.' model's 
#' estimated coefficients and frailty effects, for each unit in each time interval (represented by its mid point). 
#'
#' @param result S3 object of class 'AdPaik' containing model results.
#' @param lwd The line width of the plot. Default is 1.
#' @param xlim A numeric vector specifying the range for the x-axis (intervals). Default is min-max value of the time domain.
#' @param ylim A numeric vector specifying the range for the y-axis (intervals). Default is the range 0-1.
#' @param xlab,ylab String giving the x and y axis name. Default values are 'Time' and 'Values'.
#' @param main Title of the plot. Default title is 'Survival'.
#' @param cex Dimension of the points used for plotting the estimates. Defaults to 0.2.
#' @param cexlegend Dimension of the text used for the legend. Defaults to 0.9.
#'
#' @return The plot of the conditional survival function.
#'
#' @export
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
#' plot_survivalAdPaik(result)
#'  } 
plot_survivalAdPaik <- function(result, lwd = 1, 
                          xlim = c(min(result$TimeDomain), max(result$TimeDomain)), ylim = c(0,1),
                          xlab = "Time", ylab = "Values", main = "Conditional Survival",
                          cex = 0.2, cexlegend = 0.8){
  
  survival_df = survivalAdPaik(result)
  
  time = result$TimeDomain
  
  midpoints <- (result$TimeDomain[-1] + result$TimeDomain[-length(result$TimeDomain)]) / 2
  
  # Assuming the first column in survival_df contains the group variable
  group_variable <- survival_df[, 1]  # Extract the group variable
  group_data <- survival_df[, 2:ncol(survival_df)]  # Exclude the group variable
  
  # Order group levels using levels(factor())
  ordered_levels <- levels(factor(group_variable))  # Get ordered levels of the group variable
  set.seed(1)
  group_colors <- stats::setNames(sample(grDevices::colors(), length(ordered_levels)), ordered_levels)  # Assign colors to ordered levels
  
  n_intervals <- result$NIntervals
  
  plot(midpoints, group_data[1, ], cex = cex,
       main = main, xlab = xlab, ylab = ylab,
       xlim = xlim, ylim = ylim, col = group_colors[as.character(group_variable[1])])
  for(k in 1:(n_intervals-1))
    lines(c(midpoints[k],midpoints[k+1]), c(group_data[1,k], group_data[1,k+1]), 
          col = group_colors[as.character(group_variable[1])], lwd = lwd)
  for(i in 2:nrow(group_data)){
    for(k in 1:(n_intervals-1)){
      points(midpoints[k], group_data[i,k], cex = cex, col = group_colors[as.character(group_variable[i])])
      points(midpoints[k+1], group_data[i,k+1], cex = cex, col = group_colors[as.character(group_variable[i])])
      lines(c(midpoints[k],midpoints[k+1]), c(group_data[i,k], group_data[i,k+1]), 
            col = group_colors[as.character(group_variable[i])], lwd = lwd )
    }
  }  
  
  # Add a horizontal legend outside the plot on the right, split into two rows
  legend(
    "bottomleft",  
    legend = names(group_colors),  # Ordered group labels
    col = group_colors, 
    lwd = lwd, 
    title = "Groups", 
    cex = cexlegend,
    ncol = 3      # Split the legend into columns
  )
}




