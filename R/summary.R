#' @title
#' Summary of the 'Adapted Paik et al.'s Model'
#'
#' @description
#' Summary function for summarizing the most important information related to the dataset (number of individuals,
#' number of regressors, number of intervals, number of clusters), the model call (number of parameters) and the model
#' output (optimal log-likelihood value and AIC).
#'
#' @details
#' Among the estimated parameters, only the regressors are reported together with their standard error and confidence interval.
#'
#' @param result 'S3' class object returned by the main model call, i.e. output of the 'Adapted Paik et al.'s Model'.
#'
#' @return Model summary printed on output.
#'
#' @export
#'
#' @examples
#' # Define the variables needed for the model execution
#' formula <-
#' time_axis <-
#' categories_range_min <- c()
#' categories_range_max <- c()
#'
#' # Call the main model function
#' result <- AdPaikModel(formula, data, time_axis, categories_range_min, categories_range_max, flag_fullsd = TRUE)
#'
#' # Call the summary
#' summary.AdPaik(result)

summary.AdPaik <- function(result){
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
  string_parameters <- paste(string_parameters, "divided as (phi, betar, mu1, nu, gammak) = (", sep=",\n")
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
  
  output <- paste("Output of the 'Adapted Paik et al.'s Model'", paste4, sep="\n")
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
#' @title
#' Summary of the 'Adapted Paik et al.'s Model'
#'
#' @description
#' Summary function for summarizing the most important information related to the dataset (number of individuals,
#' number of regressors, number of intervals, number of clusters), the model call (number of parameters) and the model
#' output (optimal log-likelihood value and AIC).
#'
#' @details
#' Among the estimated parameters, only the regressors are reported together with their standard error and confidence interval.
#'
#' @param result 'S3' class object returned by the main model call, i.e. output of the 'Adapted Paik et al.'s Model'.
#'
#' @return Model summary printed on output.
#'
#' @export
#'
#' @examples
#' # Define the variables needed for the model execution
#' formula <-
#' time_axis <-
#' categories_range_min <- c()
#' categories_range_max <- c()
#'
#' # Call the main model function
#' result <- AdPaikModel(formula, data, time_axis, categories_range_min, categories_range_max, flag_fullsd = TRUE)
#'
#' # Call the summary
#' summary.AdPaik(result)

summary.StocTimeDep <- function(result){
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
  string_parameters <- paste(string_parameters, "divided as (phi, betar, sigmac, sigmab, sigmacb) = (", sep=",\n")
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
  
  output <- paste("Output of the 'Stochastic Time-Dependent Centre-Specific Frailty Model'", paste4, sep="\n")
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
#' Function for plotting on the R console a summary of the model call, given the
#' output of the same call.
#' According to the model it has been called, this function inspects the class output and then call
#' the proper summary function.
#'
#' @result Output of the model call. The result has its own class, associated to the called model.

summary <- function(result){
  if(class(result) == "AdPaik")
    summary.AdPaik(result)
  else if(class(result) == "PowPar")
    summary.PowPar(result)
  else if(class(result) == "StocTimeDep")
    summary.StocTimeDep(result)
}

#-------------------------------------------------------------------------------
