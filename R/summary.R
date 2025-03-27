#' @title Summary Method for AdPaik Objects
#'
#' @description
#' `summary` method for objects of class `"AdPaik"`. Provides a structured summary of the results from the Adapted Paik et al.'s Time-Dependent Shared Frailty Model.
#'
#' @details
#' This function extracts relevant model outputs, such as estimated parameters, standard errors,
#' and the convergence status of the algorithm. The summary is returned as an object of class `"summary.AdPaik"`,
#' which can be printed using `print()`.
#'
#' @param object An object of class `"AdPaik"`, returned by the main model function.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' `summary.AdPaik()` returns an object of class `"summary.AdPaik"`, containing structured model summary information.
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
#' summary_result <- summary(result)
#' print(summary_result)  # or simply `summary(result)`
#' }
summary.AdPaik <- function(object, ...) {
  check.result(object) # Validate input
  
  # Extract model information
  params_categories <- object$ParametersCategories
  L <- params_categories[1] # Number of intervals
  R <- params_categories[2] # Number of regressors
  
  # Format optimal parameters with standard errors
  n_params <- object$NParameters
  optimal_parameters <- sapply(1:n_params, function(p) {
    sprintf("%.4f (%.4f)", object$OptimalParameters[p], object$StandardErrorParameters[p])
  })
  
  # Extract estimated regressors
  betar <- optimal_parameters[(L + 1):(L + R)]
  
  # Convergence status
  convergence <- if (object$Status) {
    sprintf("TRUE (Converged in %d runs)", object$NRun)
  } else {
    "FALSE (No Convergence)"
  }
  
  # Construct structured summary output
  summary_list <- list(
    call = as.character(object$formula),
    cluster_variable = object$ClusterVariable,
    n_clusters = object$NClusters,
    log_likelihood = round(object$Loglikelihood, 4),
    AIC = round(object$AIC, 4),
    convergence = convergence,
    n_parameters = object$NParameters,
    parameter_breakdown = params_categories,
    n_intervals = L,
    n_regressors = R,
    regressors = stats::setNames(betar, object$Regressors)
  )
  
  class(summary_list) <- "summary.AdPaik"
  return(summary_list)
}


#' @title Print Method for summary.AdPaik Objects
#'
#' @description
#' `print` method for objects of class `"summary.AdPaik"`. Prints the formatted summary to the console.
#'
#' @param x An object of class `"summary.AdPaik"`, created by `summary(object)`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' `print.summary.AdPaik()` prints the formatted summary to the console.
#'
#' @seealso \code{\link{summary}}, \code{\link{print}}
#' @export
#' @method print summary.AdPaik
print.summary.AdPaik <- function(x, ...) {
  cat("Output of the 'Adapted Paik et al.' Model\n")
  cat("---------------------------------------------------\n")
  cat("Call: ", paste(x$call, collapse = " "), "\n")
  cat(sprintf("Cluster variable: '%s' (%d clusters)\n", x$cluster_variable, x$n_clusters))
  cat("---------------------------------------------------\n")
  cat(sprintf("Log-likelihood: %.4f\n", x$log_likelihood))
  cat(sprintf("AIC: %.4f\n", x$AIC))
  cat(sprintf("Status of the algorithm: %s\n", x$convergence))
  cat("---------------------------------------------------\n")
  cat(sprintf("Overall number of parameters: %d\n", x$n_parameters))
  cat("Parameters breakdown: ", paste(x$parameter_breakdown, collapse = ", "), "\n")
  cat(sprintf("Number of intervals: %d\n", x$n_intervals))
  cat(sprintf("Number of coefficients: %d\n", x$n_regressors))
  cat("---------------------------------------------------\n")
  cat("Estimated coefficients (with standard errors):\n")
  for (r in seq_along(x$regressors)) {
    cat(sprintf("  %s: %s\n", names(x$regressors)[r], x$regressors[r]))
  }
  cat("---------------------------------------------------\n")
}
