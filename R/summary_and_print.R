#' @title Summary Method for AdPaik Objects
#'
#' @description
#' `summary` method for objects of class `"AdPaik"`. It prepares a structured summary of the results.
#'
#' @param object An object of class `"AdPaik"`, returned by the main model function.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' An object of class `"summary.AdPaik"` containing structured model summary information.
#'
#' @export
#' @method summary AdPaik
summary.AdPaik <- function(object, ...) {
  check.result(object)
  
  params_categories <- object$ParametersCategories
  n_categories <- length(params_categories)
  L <- n_intervals <- params_categories[1]
  R <- n_regressors <- params_categories[2]
  
  n_params <- object$NParameters
  optimal_parameters <- rep(0, n_params)
  for (p in 1:n_params) {
    optimal_parameters[p] <- paste0(round(object$OptimalParameters[p], 4), " (", round(object$StandardErrorParameters[p], 4), ")")
  }
  
  betar <- optimal_parameters[(L + 1):(L + R)]
  
  convergence <- if (object$Status) {
    paste("TRUE (Convergence in", object$NRun, "runs).")
  } else {
    "FALSE (No Convergence)"
  }
  
  summary_list <- list(
    call = paste(object$formula[2], object$formula[1], object$formula[3]),
    cluster = paste0("Cluster variable '", object$ClusterVariable, "' (", object$NClusters, " clusters)."),
    logLik = round(object$Loglikelihood, 4),
    AIC = round(object$AIC, 4),
    convergence = convergence,
    parameters_info = paste0(
      "Overall number of estimated parameters ", object$NParameters,
      " divided as (phi, beta, mu1, nu, gammak) = (", paste(params_categories, collapse = ","), ")"
    ),
    n_intervals = n_intervals,
    n_regressors = n_regressors,
    regressors = stats::setNames(betar, object$Regressors)
  )
  
  class(summary_list) <- "summary.AdPaik"
  return(summary_list)
}

#-------------------------------------------------------------------------------

#' @title Print Method for summary.AdPaik Objects
#'
#' @description
#' `print` method for objects of class `"summary.AdPaik"`. Formats and prints the model summary.
#'
#' @param x An object of class `"summary.AdPaik"`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Prints a structured summary of the `AdPaik` model.
#'
#' @export
#' @method print summary.AdPaik
print.summary.AdPaik <- function(x, ...) {
  sep_line <- "-------------------------------------------------------------------------------"
  
  output <- c(
    "Output of the 'Adapted Paik et al.'s Model'",
    sep_line,
    paste("Call:", x$call),
    x$cluster,
    sep_line,
    paste("Log-likelihood:", x$logLik),
    paste("AIC:", x$AIC),
    paste("Status of the algorithm:", x$convergence),
    sep_line,
    x$parameters_info,
    paste("with: number of intervals =", x$n_intervals, ", number of regressors =", x$n_regressors, "."),
    sep_line,
    "Estimated regressors (standard error):"
  )
  
  for (r in names(x$regressors)) {
    output <- c(output, paste(r, ":", x$regressors[r]))
  }
  
  output <- c(output, sep_line)
  cat(paste(output, collapse = "\n"), "\n")
}




#-------------------------------------------------------------------------------
#' Print method for AdPaik objects
#'
#' This function prints a summary of the optimal parameters estimated in an AdPaik object.
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
  cat(sprintf("divided as (phi, beta, mu1, nu, gammak) = (%s)", paste(x$ParametersCategories, collapse = ",")))
  
  # Print estimated parameters
  cat("\nOptimal Parameters:\n")
  print(x$OptimalParameters)
  
  # Print standard errors if available
  if (!is.null(x$StandardErrorParameters)) {
    cat("\nStandard Errors:\n")
    print(x$StandardErrorParameters)
  }
  
  # Print standard errors if available
  if (!is.null(x$ParametersCI$ParamsCI_left)) {
    ci <- cbind(x$ParametersCI$ParamsCI_left, x$ParametersCI$ParamsCI_right)
    a <- (1 - 0.95) / 2
    pct <- format_perc(c(a, 1 - a), 3)
    colnames(ci) <- pct
    cat("\n95% Confidence Intervals:\n")
    print(ci)
  }
  
  # Return object invisibly
  invisible(x)
}
