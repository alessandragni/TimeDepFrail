#' @title
#' Numerical Check of First-Order Optimality Conditions
#'
#' @description
#' Function for computing a finite-difference approximation of the gradient of the
#' log-likelihood function at a given parameter vector (typically the optimal parameters
#' returned by the multi-dimensional Powell/Brent optimization), together with two summary
#' diagnostics: the Euclidean norm and the maximum absolute component of the gradient,
#' restricted to parameters not sitting close to one of their declared range boundaries.
#'
#' @details
#' Since the log-likelihood function is maximized using a derivative-free, coordinate-wise
#' optimization strategy, no analytical gradient is available during optimization. This
#' function provides a numerical check, computed once at the end of the optimization
#' procedure, of the first-order optimality conditions: at a genuine interior stationary
#' point, all partial derivatives of the log-likelihood should be close to zero, as opposed
#' to a flat region of the likelihood surface where the coordinate-wise optimizer could
#' stall.
#'
#' Each partial derivative is approximated through a centered finite-difference scheme. The
#' step used is always symmetric (identical forward and backward), shrunk if necessary to
#' 'min(h_grad, distance to the nearer declared range bound)': this avoids the numerical
#' instability that an asymmetric step (a full step on one side, a boundary-clamped, much
#' smaller step on the other) would introduce. If the parameter sits at (or past) a declared
#' boundary, leaving no room for any step, the corresponding entry of 'Gradient' is 'NA'. If
#' perturbing a parameter pushes the log-likelihood itself to a non-finite value on one side
#' (which can happen even strictly inside the declared range, since some of the gamma-ratio
#' terms of the log-likelihood become steep for near-degenerate dispersion parameters), a
#' one-sided difference on the finite side is used instead, or 'NA' if both sides are
#' non-finite.
#'
#' Parameters whose optimum lies within 'boundary_factor * h_grad' of either declared bound
#' are reported separately in 'BoundaryAdjacent', rather than folded into 'GradientNorm'/
#' 'GradientMaxAbs'. Such parameters (e.g. a dispersion parameter shrinking towards a lower
#' range bound of virtually zero) can show a large, genuine, one-sided gradient at a boundary
#' solution - a known, different phenomenon from an interior optimizer stall - so keeping
#' them separate avoids a boundary solution being mistaken for evidence against convergence.
#'
#' @param optimal_params Numerical vector of parameters at which the gradient is evaluated.
#' Its length (i.e. number of parameters) is equal to \eqn{n_p}.
#' @param params_range_min Numerical vector of length equal to \eqn{n_p}, containing the minimum range of each parameter.
#' @param params_range_max Numerical vector of length equal to \eqn{n_p}, containing the maximum range of each parameter.
#' @param dataset Dataset containing the value of the regressors for all individuals in the study.
#' @param centre vector containing the group membership of each individual and that induces the clustering subdivision.
#' @param time_axis Temporal domain.
#' Its number of intervals corresponds to the length of the time-domain minus 1
#' @param dropout_matrix Binary matrix of dimension (n_individuals, n_intervals).
#' The sum of the elements of each row must be (1), if the associated individual failed in a precise interval, and (0) if the individual
#' did not fail in the time-axis.
#' Therefore, if an individual failed in the time-domain, the interval in which he failed will have value (1) and the others (0).
#' @param e_matrix Matrix of dimension (n_individuals, n_intervals) where each element contains the resolution of the temporal
#' integral for that individual in that interval, thorugh the 'e_time_fun' function.
#' @param h_grad Discretization step for the numerical approximation of the first derivative (gradient) of the log-likelihood function.
#' @param boundary_factor Multiple of 'h_grad' used as the distance-to-boundary threshold for flagging a parameter as
#' 'boundary-adjacent' (see Details). Defaults to 10.
#'
#' @return A list composed of:
#' - Gradient: numerical vector of length equal to the number of model parameters, containing the finite-difference approximation
#' of each partial derivative of the log-likelihood function, evaluated at 'optimal_params' ('NA' where it could not be estimated).
#' - BoundaryAdjacent: integer vector of indices (possibly empty) of parameters whose optimum lies within 'boundary_factor * h_grad'
#' of either declared range bound; excluded from 'GradientNorm'/'GradientMaxAbs'.
#' - GradientNorm: Euclidean (L2) norm of 'Gradient', restricted to non-boundary-adjacent, non-'NA' components.
#' - GradientMaxAbs: maximum absolute component of 'Gradient', restricted to non-boundary-adjacent, non-'NA' components.
#'
#' @keywords internal

gradient_check <- function(optimal_params, params_range_min, params_range_max,
                           dataset, centre, time_axis, dropout_matrix, e_matrix, h_grad,
                           boundary_factor = 10){

  # Extract information from input variables
  n_params <- length(optimal_params)

  # Initialize gradient vector
  gradient <- rep(NA_real_, n_params)

  # Log-likelihood at the provided parameter vector: needed whenever a one-sided
  # difference is used
  ll_eval <- ll_AdPaik_eval(optimal_params, dataset, centre, time_axis, dropout_matrix, e_matrix)

  # Distance from each parameter to its nearer declared range boundary
  dist_to_bound <- pmin(optimal_params - params_range_min, params_range_max - optimal_params)

  # Compute gradient components
  for(p in 1:n_params){
    value <- optimal_params[p]

    # Symmetric step: shrink to whatever room is available on the tighter side, so
    # both sides of the finite difference always use the identical step size (this
    # avoids the instability of blending two very different step sizes when a
    # parameter sits close to one of its declared bounds).
    h_p <- min(h_grad, value - params_range_min[p], params_range_max[p] - value)

    if(h_p <= 0){
      # Parameter sits at (or past) a declared boundary: no room for any step.
      next
    }

    params_plus  <- optimal_params
    params_minus <- optimal_params
    params_plus[p]  <- value + h_p
    params_minus[p] <- value - h_p

    ll_plus  <- ll_AdPaik_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_minus <- ll_AdPaik_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)

    if(is.finite(ll_plus) && is.finite(ll_minus)){
      gradient[p] <- (ll_plus - ll_minus) / (2 * h_p)
    }
    else if(is.finite(ll_plus)){
      # Only the forward step is usable: one-sided forward difference
      gradient[p] <- (ll_plus - ll_eval) / h_p
    }
    else if(is.finite(ll_minus)){
      # Only the backward step is usable: one-sided backward difference
      gradient[p] <- (ll_eval - ll_minus) / h_p
    }
    # else: log-likelihood non-finite on both sides, leave gradient[p] as NA
  }

  # Parameters whose optimum lies close to a declared range boundary: reported
  # separately, excluded from the interior gradient diagnostics below
  is_boundary <- dist_to_bound < (boundary_factor * h_grad)
  boundary_adjacent <- which(is_boundary)

  interior_gradient <- gradient[!is_boundary]
  interior_gradient <- interior_gradient[is.finite(interior_gradient)]

  gradient_norm <- if(length(interior_gradient) == 0) NA_real_ else sqrt(sum(interior_gradient^2))
  gradient_max_abs <- if(length(interior_gradient) == 0) NA_real_ else max(abs(interior_gradient))

  return(list("Gradient" = gradient,
             "BoundaryAdjacent" = boundary_adjacent,
             "GradientNorm" = gradient_norm,
             "GradientMaxAbs" = gradient_max_abs))
}
