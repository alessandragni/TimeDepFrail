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
#' point, all partial derivatives of the log-likelihood should be close to zero.
#'
#' With at least 'h_grad' of room on both sides of its declared range, a parameter's
#' derivative is approximated by the standard centered finite difference with step 'h_grad'.
#' If the optimum lies within 'h_grad' of a bound, a symmetric step shrunk to match that
#' distance is deliberately *not* used: both evaluation points would then straddle the
#' (possibly boundary) optimum at a resolution finer than the Brent search that reported it
#' can distinguish, so the result is dominated by curvature noise rather than a real slope
#' (verified empirically to give inflated, inconsistently-signed values). Instead, a
#' one-sided difference with the *full*, unshrunk 'h_grad' step is used, in whichever
#' direction has that much room - also the statistically meaningful quantity at a boundary
#' (the one-directional score). Only if the whole declared range is narrower than 'h_grad'
#' does the function fall back to a symmetric step using whatever small room is available.
#' 'NA' if no room exists on either side, or if perturbing the parameter makes the
#' log-likelihood non-finite on the evaluated side(s).
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
    lower_room <- value - params_range_min[p]
    upper_room <- params_range_max[p] - value

    if(lower_room >= h_grad && upper_room >= h_grad){
      # Standard interior case: full symmetric centered difference.
      params_plus  <- optimal_params
      params_minus <- optimal_params
      params_plus[p]  <- value + h_grad
      params_minus[p] <- value - h_grad

      ll_plus  <- ll_AdPaik_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
      ll_minus <- ll_AdPaik_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)

      if(is.finite(ll_plus) && is.finite(ll_minus)){
        gradient[p] <- (ll_plus - ll_minus) / (2 * h_grad)
      }
      else if(is.finite(ll_plus)){
        gradient[p] <- (ll_plus - ll_eval) / h_grad
      }
      else if(is.finite(ll_minus)){
        gradient[p] <- (ll_eval - ll_minus) / h_grad
      }
    }
    else if(upper_room >= h_grad || lower_room >= h_grad){
      # Boundary-adjacent: a symmetric step would have to be shrunk well below
      # h_grad on the near side, which produces a numerically meaningless
      # result at this resolution (see Details). Use a one-sided difference
      # with the full, unshrunk h_grad step in the feasible direction instead.
      params_step <- optimal_params
      if(upper_room >= h_grad){
        params_step[p] <- value + h_grad
        ll_step <- ll_AdPaik_eval(params_step, dataset, centre, time_axis, dropout_matrix, e_matrix)
        if(is.finite(ll_step)) gradient[p] <- (ll_step - ll_eval) / h_grad
      }
      else{
        params_step[p] <- value - h_grad
        ll_step <- ll_AdPaik_eval(params_step, dataset, centre, time_axis, dropout_matrix, e_matrix)
        if(is.finite(ll_step)) gradient[p] <- (ll_eval - ll_step) / h_grad
      }
    }
    else if(lower_room > 0 && upper_room > 0){
      # Degenerate range, narrower than h_grad on both sides: fall back to a
      # symmetric step using whatever (small) room is available on both sides
      # equally. Reduced accuracy, but keeps this edge case well-defined.
      h_p <- min(lower_room, upper_room)
      params_plus  <- optimal_params
      params_minus <- optimal_params
      params_plus[p]  <- value + h_p
      params_minus[p] <- value - h_p

      ll_plus  <- ll_AdPaik_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
      ll_minus <- ll_AdPaik_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)

      if(is.finite(ll_plus) && is.finite(ll_minus)){
        gradient[p] <- (ll_plus - ll_minus) / (2 * h_p)
      }
    }
    # else: no room on either side (parameter sits at, or past, a declared
    # boundary) -- gradient[p] stays NA.
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
