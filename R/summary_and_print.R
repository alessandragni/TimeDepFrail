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

  gradient_norm <- round(object$GradientCheck$GradientNorm, 6)
  gradient_max_abs <- round(object$GradientCheck$GradientMaxAbs, 6)
  n_boundary <- length(object$GradientCheck$BoundaryAdjacent)
  boundary_note <- if (n_boundary == 0) {
    "none"
  } else {
    paste0(n_boundary, " parameter(s) at index/indices ",
           paste(object$GradientCheck$BoundaryAdjacent, collapse = ", "),
           " excluded (boundary-adjacent)")
  }
  gradient_check <- paste0("||grad||_2 = ", gradient_norm, ", max|grad| = ", gradient_max_abs,
                          "\n[", boundary_note, "]")

  # Corrected (full-Hessian, cross-derivative-aware) regressor standard errors,
  # only available if 'full_hessian_se = TRUE' was passed to AdPaikModel()
  regressors_corrected <- NULL
  if(!is.null(object$HessianCheck)){
    se_full_betar <- object$HessianCheck$SE_full[(L + 1):(L + R)]
    if(!anyNA(se_full_betar)){
      regressors_corrected <- stats::setNames(
        paste0(round(se_full_betar, 4), " (was ", round(object$StandardErrorParameters[(L + 1):(L + R)], 4), ")"),
        object$Regressors)
    }
  }

  summary_list <- list(
    call = paste(object$formula[2], object$formula[1], object$formula[3]),
    cluster = paste0("Cluster variable '", object$ClusterVariable, "' (", object$NClusters, " clusters)."),
    logLik = round(object$Loglikelihood, 4),
    AIC = round(object$AIC, 4),
    convergence = convergence,
    gradient_check = gradient_check,
    parameters_info = paste0(
      "Overall number of estimated parameters ", object$NParameters,
      " divided as (phi, beta, mu1, nu, gammak) = (", paste(params_categories, collapse = ","), ")"
    ),
    n_intervals = n_intervals,
    n_regressors = n_regressors,
    regressors = stats::setNames(betar, object$Regressors),
    regressors_corrected = regressors_corrected
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
    paste("First-order optimality check (finite-difference gradient):", x$gradient_check),
    sep_line,
    x$parameters_info,
    paste("with: number of intervals =", x$n_intervals, ", number of regressors =", x$n_regressors, "."),
    sep_line,
    "Estimated regressors (standard error):"
  )
  
  for (r in names(x$regressors)) {
    output <- c(output, paste(r, ":", x$regressors[r]))
  }

  if (!is.null(x$regressors_corrected)) {
    output <- c(output, sep_line, "Corrected (full-Hessian, cross-derivative-aware) standard error:")
    for (r in names(x$regressors_corrected)) {
      output <- c(output, paste(r, ":", x$regressors_corrected[r]))
    }
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

  # Print first-order optimality check (finite-difference gradient at the optimum)
  if (!is.null(x$GradientCheck)) {
    cat("\nFirst-order optimality check (finite-difference gradient at the optimum):\n")
    cat(sprintf("  Interior parameters: ||grad||_2 = %s, max|grad| = %s\n",
                signif(x$GradientCheck$GradientNorm, 6),
                signif(x$GradientCheck$GradientMaxAbs, 6)))
    boundary_idx <- x$GradientCheck$BoundaryAdjacent
    if (length(boundary_idx) > 0) {
      cat("  Boundary-adjacent parameters (excluded above; large one-sided gradient expected):\n")
      cat(sprintf("    index %d: grad = %s\n", boundary_idx, signif(x$GradientCheck$Gradient[boundary_idx], 6)), sep = "")
    }
  }

  # Print full-Hessian, correlation-aware standard errors, if requested via
  # 'full_hessian_se = TRUE' in the AdPaikModel() call
  if (!is.null(x$HessianCheck)) {
    cat("\nFull-Hessian (cross-derivative-aware) standard errors:\n")
    cat(sprintf("  Interior covariance matrix positive definite: %s\n", x$HessianCheck$IsPositiveDefinite))
    cat("  SE_full (interior parameters; NA = boundary-adjacent):\n")
    print(x$HessianCheck$SE_full)
    cat("  SE_full / StandardErrorParameters ratio:\n")
    print(x$HessianCheck$SE_ratio)
    if (!is.na(x$HessianCheck$CorrelationMu1Nu)) {
      cat(sprintf("  Correlation(mu1, nu) = %s\n", signif(x$HessianCheck$CorrelationMu1Nu, 4)))
    } else {
      cat("  Correlation(mu1, nu): not available (raw 2x2 Hessian sub-block is not negative definite)\n")
    }
  }

  # Return object invisibly
  invisible(x)
}
