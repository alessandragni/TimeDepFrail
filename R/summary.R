#' @title
#' Summary of the Adapted Paik et al.'s Time-Dependent Shared Frailty Model
#'
#' @description
#' This function provides a comprehensive summary of the results from the Adapted Paik et al.'s Time-Dependent Shared Frailty Model. 
#' It includes key information about the dataset (e.g., number of individuals, regressors, intervals, and clusters), model parameters, 
#' and output (log-likelihood, AIC). The summary also lists the estimated regressors along with their standard errors
#' 
#' @details
#' The function reports the estimated regressors, their standard errors, and confidence intervals (if available). 
#'
#' @method summary AdPaik
#' 
#' @param result 'S3' class object returned by the main model call, i.e. output of the 'Adapted Paik et al.'s Model'.
#' 
#' @return Model summary printed on output.
#'
#' @export
#'
#' @examples
#' # Define the variables needed for the model execution
#' 
#' data(data_dropout)
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
#' # Call the summary
#' summary(result)
#' }

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
#' @title Summary for Time-Dependent Frailty Models
#' 
#' @description
#' This function displays a summary of the model output based on the class of the result object. 
#' It delegates to the appropriate summary method according to the class of the result.
#'
#' @param result An object containing the output of the model call. 
#' The class of this object determines which summary method is invoked.
#' 
#' @return A summary of the model output printed to the console.
#' 
#' @export
summary <- function(result){
  if(inherits(result, "AdPaik"))
    summary.AdPaik(result)
}

#-------------------------------------------------------------------------------
