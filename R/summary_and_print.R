#' @title Summary Method for AdPaik Objects
#'
#' @description
#' `summary` method for objects of class `"AdPaik"`. Provides a structured summary of the results from the Adapted Paik et al.'s Time-Dependent Shared Frailty Model.
#'
#' @details
#' This function extracts relevant model outputs, such as estimated parameters, standard errors,
#' and the convergence status of the algorithm. 
#'
#' @param object An object of class `"AdPaik"`, returned by the main model function.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' `summary.AdPaik()` returns a structured model summary information.
#'
#' @seealso \code{\link{summary}}, \code{\link{print}}
#' @export
#' @method summary AdPaik
#'
#' @examples
#' \dontrun{
#' data(data_dropout)
#' eps_paik <- 1e-10
#' categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
#' categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
#' time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
#' formula <- time_to_event ~ Gender + CFUP + cluster(group)
#'
#' # Fit the model
#' result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#'
#' # Generate and print summary
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
  
  string_parameters <- paste("Overall number of estimated parameters ", result$NParameters)
  string_parameters <- paste(string_parameters, "divided as (phi, beta, mu1, nu, gammak) = (", sep=",\n")
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
  cat("\n")
}




#-------------------------------------------------------------------------------
#' Print method for AdPaik objects
#'
#' This function prints a summary of an AdPaik object.
#'
#' @param x An object of class `AdPaik`.
#' @param ... Additional arguments (not used).
#'
#' @return Prints a summary of the `AdPaik` object and returns it invisibly.
#' @export
print.AdPaik <- function(x, ...) {
  # Check if x is an AdPaik object
  if (!inherits(x, "AdPaik")) {
    stop("The input object is not of class 'AdPaik'.")
  }
  
  cat("\n'Adapted Paik et al.' Model\n")
  cat(rep("-", 30), "\n", sep = "")
  
  # Print basic model structure
  # Fix: Use deparse() to handle formula safely
  if (!is.null(x$formula)) {
    cat("Formula: ", paste(deparse(x$formula), collapse = " "), "\n")
  }
  
  cat("Number of Intervals: ", x$NIntervals, "\n")
  cat("Number of Regressors:", x$NRegressors, "\n")
  
  # Print regressors if available
  if (!is.null(x$Regressors) && length(x$Regressors) > 0) {
    cat("Regressors: ", paste(x$Regressors, collapse = ", "), "\n")
  }
  
  # Fix incorrect usage of cat() with formatted printing
  cat(sprintf("Overall number of parameters: %d\n", x$NParameters))
  
  # Print estimated parameters
  cat("\nEstimated Parameters:\n")
  print(x$OptimalParameters)
  
  # Print standard errors if available
  if (!is.null(x$StandardErrorParameters)) {
    cat("\nStandard Errors:\n")
    print(x$StandardErrorParameters)
  }
  
  # Return object invisibly
  invisible(x)
}
