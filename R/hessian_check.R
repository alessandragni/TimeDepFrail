#' @title
#' Full Numerical Hessian and Correlation-Aware Standard Errors
#'
#' @description
#' Function for computing the full (not diagonal-only) finite-difference Hessian of the
#' log-likelihood function at a given parameter vector (typically the optimal parameters
#' returned by the multi-dimensional Powell/Brent optimization), together with the
#' resulting correlation-aware standard errors, for parameters that are not boundary-adjacent.
#'
#' @details
#' `params_se()` (`R/params_se_CI.R`) only approximates the diagonal of the Hessian, for
#' computational efficiency during repeated calls; this function computes the full matrix,
#' but only once, at a fixed parameter vector (typically the already-reported optimum), so the
#' additional cost (\eqn{O(n_p^2)} extra log-likelihood evaluations instead of \eqn{O(n_p)}) is
#' paid a single time.
#'
#' Each parameter uses a symmetric, boundary-aware step `h <- min(h_dd, value - range_min,
#' range_max - value)`, matching the step-size logic already adopted in `gradient_check()`
#' (`R/gradient_check.R`) to avoid the instability of an asymmetric clamped step.
#'
#' Since a parameter sitting at/near a declared range boundary can make the full Hessian
#' ill-conditioned or non-negative-definite overall (a near-zero-curvature boundary direction
#' couples into every other row/column once the matrix is inverted, unlike the diagonal-only
#' approximation, which is immune to other parameters' degeneracy by construction), the
#' covariance matrix and its standard errors are computed only from the submatrix of
#' non-boundary-adjacent ("interior") parameters. Boundary-adjacent parameters (see
#' `gradient_check()`) keep their diagonal-only standard error (from `params_se()`) unchanged;
#' comparing the two is only meaningful for interior parameters.
#'
#' @param optimal_params Numerical vector of parameters at which the Hessian is evaluated.
#' Its length (i.e. number of parameters) is equal to \eqn{n_p}.
#' @param params_range_min Numerical vector of length equal to \eqn{n_p}, containing the minimum range of each parameter.
#' @param params_range_max Numerical vector of length equal to \eqn{n_p}, containing the maximum range of each parameter.
#' @param dataset Dataset containing the value of the regressors for all individuals in the study.
#' @param centre vector containing the group membership of each individual and that induces the clustering subdivision.
#' @param time_axis Temporal domain.
#' @param dropout_matrix Binary matrix of dimension (n_individuals, n_intervals).
#' @param e_matrix Matrix of dimension (n_individuals, n_intervals), see `params_se()`.
#' @param h_dd Discretization step for the numerical approximation of the second derivative of the log-likelihood function.
#' @param boundary_adjacent Integer vector of indices of boundary-adjacent parameters (as returned
#' by `gradient_check()$BoundaryAdjacent`), excluded from the interior covariance/SE computation.
#' @param se_diag Numerical vector of the existing diagonal-only standard errors (from `params_se()`),
#' used unchanged for boundary-adjacent parameters and for comparison for interior ones.
#' @param mu1_index,nu_index Position of `mu1` and `nu` in the parameter vector. Their raw 2x2
#' Hessian sub-block is reported separately (see Details), since Section 3.3's `mu1`/`nu` algebraic
#' coupling makes this pair of particular interest regardless of whether it ends up interior or
#' boundary-adjacent, where a full-matrix inversion may not be available.
#'
#' @return A list composed of:
#' - Hessian: full \eqn{n_p \times n_p} numerical Hessian matrix (`NA` entries where a boundary
#' left no room for the required step).
#' - InteriorIndices: integer vector of the non-boundary-adjacent parameter indices actually used
#' for the covariance/SE computation below.
#' - SE_full: numerical vector of length \eqn{n_p}: correlation-aware standard errors for interior
#' parameters, `NA` for boundary-adjacent ones.
#' - SE_ratio: `SE_full / se_diag` for interior parameters (`NA` for boundary-adjacent), summarizing
#' how much the full Hessian changes each interior parameter's standard error relative to the
#' diagonal-only approximation.
#' - CorrelationMu1Nu: the (`mu1`,`nu`) correlation implied by the raw 2x2 `mu1`/`nu` Hessian
#' sub-block (`NA` if that sub-block is not positive definite), reported independently of the
#' interior/boundary-adjacent split, since Section 3.3's `mu1`/`nu` coupling makes this pair of
#' particular interest regardless of which side of that split it falls on.
#' - Mu1NuBlockPositiveDefinite: logical, whether the raw 2x2 `mu1`/`nu` Hessian sub-block is
#' positive definite (i.e. `CorrelationMu1Nu` is a valid correlation and not `NA`).
#' - IsPositiveDefinite: logical, whether `-Hessian[InteriorIndices, InteriorIndices]` is positive
#' definite (i.e. the interior covariance matrix is a valid one).
#'
#' @keywords internal

hessian_check <- function(optimal_params, params_range_min, params_range_max,
                          dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd,
                          boundary_adjacent, se_diag, mu1_index, nu_index){

  n_params <- length(optimal_params)

  # Symmetric, boundary-aware step per parameter (same logic as gradient_check())
  h <- pmin(h_dd, optimal_params - params_range_min, params_range_max - optimal_params)

  ll0 <- ll_AdPaik_eval(optimal_params, dataset, centre, time_axis, dropout_matrix, e_matrix)

  # Single-parameter perturbations, reused for both diagonal and off-diagonal terms
  ll_plus <- ll_minus <- rep(NA_real_, n_params)
  for(p in 1:n_params){
    if(h[p] <= 0) next
    params_plus <- params_minus <- optimal_params
    params_plus[p] <- optimal_params[p] + h[p]
    params_minus[p] <- optimal_params[p] - h[p]
    ll_plus[p] <- ll_AdPaik_eval(params_plus, dataset, centre, time_axis, dropout_matrix, e_matrix)
    ll_minus[p] <- ll_AdPaik_eval(params_minus, dataset, centre, time_axis, dropout_matrix, e_matrix)
  }

  H <- matrix(NA_real_, n_params, n_params)

  # Diagonal terms: centered second difference
  for(p in 1:n_params){
    if(h[p] <= 0) next
    if(is.finite(ll_plus[p]) && is.finite(ll_minus[p]))
      H[p,p] <- (ll_plus[p] - 2*ll0 + ll_minus[p]) / h[p]^2
  }

  # Off-diagonal terms: 4-point mixed partial finite difference
  for(p in 1:(n_params - 1)){
    if(h[p] <= 0) next
    for(q in (p+1):n_params){
      if(h[q] <= 0) next

      params_pp <- params_pm <- params_mp <- params_mm <- optimal_params
      params_pp[p] <- optimal_params[p] + h[p]; params_pp[q] <- optimal_params[q] + h[q]
      params_pm[p] <- optimal_params[p] + h[p]; params_pm[q] <- optimal_params[q] - h[q]
      params_mp[p] <- optimal_params[p] - h[p]; params_mp[q] <- optimal_params[q] + h[q]
      params_mm[p] <- optimal_params[p] - h[p]; params_mm[q] <- optimal_params[q] - h[q]

      ll_pp <- ll_AdPaik_eval(params_pp, dataset, centre, time_axis, dropout_matrix, e_matrix)
      ll_pm <- ll_AdPaik_eval(params_pm, dataset, centre, time_axis, dropout_matrix, e_matrix)
      ll_mp <- ll_AdPaik_eval(params_mp, dataset, centre, time_axis, dropout_matrix, e_matrix)
      ll_mm <- ll_AdPaik_eval(params_mm, dataset, centre, time_axis, dropout_matrix, e_matrix)

      if(all(is.finite(c(ll_pp, ll_pm, ll_mp, ll_mm)))){
        H[p,q] <- H[q,p] <- (ll_pp - ll_pm - ll_mp + ll_mm) / (4 * h[p] * h[q])
      }
    }
  }

  # Interior (non-boundary-adjacent) parameters: covariance/SE only computed here
  interior_indices <- setdiff(seq_len(n_params), boundary_adjacent)
  H_interior <- H[interior_indices, interior_indices, drop = FALSE]

  se_full <- rep(NA_real_, n_params)
  is_pd <- FALSE
  Sigma_interior <- NULL

  if(length(interior_indices) > 0 && !anyNA(H_interior)){
    eig <- eigen(-H_interior, symmetric = TRUE, only.values = TRUE)$values
    is_pd <- all(eig > 0)
    if(is_pd){
      Sigma_interior <- solve(-H_interior)
      se_full[interior_indices] <- sqrt(diag(Sigma_interior))
    }
  }

  se_ratio <- rep(NA_real_, n_params)
  valid <- !is.na(se_full) & !is.na(se_diag) & se_diag > 0
  se_ratio[valid] <- se_full[valid] / se_diag[valid]

  # mu1/nu correlation: reported from the raw 2x2 Hessian sub-block directly, so it
  # is available regardless of whether the pair ended up interior or
  # boundary-adjacent (where a full-matrix inversion may not exist/be valid).
  correlation_mu1_nu <- NA_real_
  mu1_nu_block_pd <- FALSE
  if(!is.na(H[mu1_index, mu1_index]) && !is.na(H[nu_index, nu_index]) &&
     !is.na(H[mu1_index, nu_index])){
    H_mu1_nu <- matrix(c(H[mu1_index, mu1_index], H[mu1_index, nu_index],
                        H[nu_index, mu1_index], H[nu_index, nu_index]), 2, 2)
    eig_mu1_nu <- eigen(-H_mu1_nu, symmetric = TRUE, only.values = TRUE)$values
    mu1_nu_block_pd <- all(eig_mu1_nu > 0)
    if(mu1_nu_block_pd){
      Sigma_mu1_nu <- solve(-H_mu1_nu)
      correlation_mu1_nu <- Sigma_mu1_nu[1,2] / sqrt(Sigma_mu1_nu[1,1] * Sigma_mu1_nu[2,2])
    }
  }

  return(list("Hessian" = H,
             "InteriorIndices" = interior_indices,
             "SE_full" = se_full,
             "SE_ratio" = se_ratio,
             "IsPositiveDefinite" = is_pd,
             "CorrelationMu1Nu" = correlation_mu1_nu,
             "Mu1NuBlockPositiveDefinite" = mu1_nu_block_pd))
}
