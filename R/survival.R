#' @title
#' Compute Survival Function
#'
#' @description
#' Computes the survival function based on the 'Adapted Paik et al.' model's 
#' estimated coefficients and frailty effects.
#'
#' @param result S3 object of class 'AdPaik' containing model results.
#' @param data Data frame containing the dataset with covariates used in the model.
#'
#' @return A dataset where each row corresponds to the survival function values 
#' over the time intervals for each individual in the dataset.
#'
#' @export
survival <- function(result, data) {
  # Check for valid input
  if (class(result) != "AdPaik") stop("'result' must be of class 'AdPaik'.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  
  # Extract beta coefficients and formula
  beta <- coef(result)$beta  # Extract coefficients
  # names(coef(result)$beta) # Extract names
  full_formula <- result$formula  # Extract the full formula
  
  # Extract covariates excluding cluster terms
  terms_object <- terms(full_formula)
  covariates <- attr(terms_object, "term.labels")
  covariates_nocluster <- covariates[!grepl("^cluster\\(", covariates)]
  
  # Construct a design matrix for covariates, accounting for factors/dummies
  covariate_data <- data[, covariates_nocluster, drop = FALSE]
  design_matrix <- model.matrix(~ . - 1, data = covariate_data)  # No intercept
  
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
  exp_phi <- exp(c(coef(result)$phi))  # Baseline hazard parameters
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
  
  # Return survival function
  return(survival)
  
}

# survival_df = survival(result, data_dropout)
# 
# library(ggplot2)
# library(reshape2)
# 
# library(tidyr, dplyr)
# library(ggplot2)
# 
# plot(result$TimeDomain, survival_df[1,], type='l', ylim = c(0,1))
# 
# 
# # Predefine colors for groups
# unique_groups <- unique(survival_df$group)
# group_colors <- setNames(rainbow(length(unique_groups)), unique_groups)
# 
# # Plot the first row as a base
# plot(result$TimeDomain, c(1, survival_df[1, -1]), type = "l",
#      ylim = c(0, 1), xlab = "Time", ylab = "Survival Probability",
#      col = group_colors[survival_df$group[1]], lwd = 2)
# 
# # Loop through remaining rows
# for (i in 2:nrow(survival_df)) {
#   lines(result$TimeDomain, c(1, survival_df[i, -1]),
#         col = group_colors[survival_df$group[i]], lwd = 2)
# }
# 
# # Add legend
# legend("bottom", legend = unique_groups, col = group_colors, lty = 1, lwd = 2, title = "Groups")




