# vcov.AdPaik <- function(optimal_params, params_range_min, params_range_max,
#                         dataset, centre, time_axis, dropout_matrix, e_matrix, h_dd) {
#   
#   n_params <- length(optimal_params)
#   H <- matrix(0, n_params, n_params)  # Hessian matrix initialization
#   
#   # Compute Hessian matrix
#   for(i in 1:n_params){
#     for(j in 1:n_params){
#       
#       # Original values
#       theta_ij <- optimal_params
#       
#       # Perturb parameters
#       theta_ij[i] <- optimal_params[i] + h_dd
#       theta_ij[j] <- optimal_params[j] + h_dd
#       ll_ij <- ll_AdPaik_eval(theta_ij, dataset, centre, time_axis, dropout_matrix, e_matrix)
#       
#       theta_i <- optimal_params
#       theta_i[i] <- optimal_params[i] + h_dd
#       ll_i <- ll_AdPaik_eval(theta_i, dataset, centre, time_axis, dropout_matrix, e_matrix)
#       
#       theta_j <- optimal_params
#       theta_j[j] <- optimal_params[j] + h_dd
#       ll_j <- ll_AdPaik_eval(theta_j, dataset, centre, time_axis, dropout_matrix, e_matrix)
#       
#       ll_orig <- ll_AdPaik_eval(optimal_params, dataset, centre, time_axis, dropout_matrix, e_matrix)
#       
#       # Compute second-order derivative approximation
#       H[i, j] <- (ll_ij - ll_i - ll_j + ll_orig) / (h_dd * h_dd)
#     }
#   }
#   
#   # Compute Fisher Information Matrix
#   I_fisher <- -H  
#   
#   # Compute variance-covariance matrix as the inverse of Fisher Information
#   if (det(I_fisher) == 0) {
#     warning("Fisher Information Matrix is singular; using pseudo-inverse.")
#     vcov_matrix <- MASS::ginv(I_fisher)  # Use pseudo-inverse for singular matrix
#   } else {
#     vcov_matrix <- solve(I_fisher)  # Inverse of Hessian
#   }
#   
#   return(vcov_matrix)
# }
