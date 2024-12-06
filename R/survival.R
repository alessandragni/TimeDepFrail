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
#' @return A matrix where each row corresponds to the survival function values 
#' over the time intervals for each individual in the dataset.
#'
#' @export
survival <- function(result, data) {
  # Check for valid input
  if (class(result) != "AdPaik") stop("'result' must be of class 'AdPaik'.")
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  
  # Extract beta coefficients and covariates
  beta <- coef(result)$beta  # Extract coefficients
  full_formula <- result$formula  # Extract the full formula
  
  # Extract covariates excluding cluster terms
  terms_object <- terms(full_formula)
  covariates <- attr(terms_object, "term.labels")
  covariates_cleaned <- covariates[!grepl("^cluster\\(", covariates)]
  
  # Construct a design matrix for covariates, accounting for factors/dummies
  covariate_data <- data[, covariates_cleaned, drop = FALSE]
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
  exp_linear_predictor <- exp(linear_predictor)
  
  # Compute the baseline hazard step-function
  phi <- exp(coef(result)$phi)  # Baseline hazard parameters
  posterior_frailty <- result$PosteriorFrailtyEstimates$Z  # Frailty estimates
  weighted_frailty <- t(t(posterior_frailty) * phi)  # Scale by hazard params
  
  # Compute cumulative hazard
  time_diffs <- diff(result$TimeDomain)  # Differences in the time domain
  cumulative_hazard <- t(apply(apply(weighted_frailty, 2, cumsum), 1, cumsum)) * time_diffs
  
  df1 = data.frame('group'=data[[result$ClusterVariable]], 'exp_linear_predictor'=exp_linear_predictor)
  df2 = data.frame('group'=levels(factor(data[[result$ClusterVariable]])), 'cumhaz'=cumulative_hazard)
  
  result_df <- merge(df1, df2, by = "group", all.x = TRUE)
  
  data.frame('group'=result_df$group,
             t(apply(- result_df$exp_linear_predictor * result_df[,3:ncol(result_df)], 1, function(row) exp(-row))))
  # Compute the survival function
  survival_function 
  
  # Return survival function
  return(survival_function)
  
}






# survival = function(result, data){
#   coef(result)$beta
#   covariates = data[i take only the covariates in coef result]
#   
#   exp(coef(result)$beta %*% covariates) 
#   
#   MyVector = exp(coef(result)$phi)
#   MyMatrix = result$PosteriorFrailtyEstimates$Z
#   MyNewMatrix = t(t(MyMatrix) * MyVector)
#   
#   integral = t(apply(apply(MyNewMatrix, 2, cumsum), 1, cumsum)) * diff(result$TimeDomain) 
#   
#   exp(- integral)
#   
#   # Extract the formula from result
#   full_formula <- result$formula
#   
#   # Get the terms from the formula
#   terms_object <- terms(full_formula)
#   
#   # Extract covariate labels, excluding "cluster"
#   covariates <- attr(terms_object, "term.labels")
#   covariates_cleaned <- covariates[!grepl("^cluster\\(", covariates)]
#   
#   # Display the cleaned covariates
#   data_dropout[covariates_cleaned]
#   
#   result$ClusterVariable
# }# 