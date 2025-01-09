
devtools::install_github("alessandragni/TimeDepFrail")

library(TimeDepFrail)
data(data_dropout)
head(data_dropout)

#-------------------------------------------------------------------------------
# MODEL APPLICATION
#-------------------------------------------------------------------------------

## ADAPTED PAIK ET AL MODEL
# Model call
eps_paik <- 1e-10
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)

time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)

formula <- time_to_event ~ Gender + CFUP + cluster(group)

result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)
summary.AdPaik(result)
summary(result)

L <- result$NIntervals
R <- result$NRegressors
result$OptimalParameters[(L+1):(L+R)]
result$StandardErrorParameters[(L+1):(L+R)]

result$OptimalParameters[L+R+2]
result$StandardErrorParameters[L+R+2]

# Plot baseline hazard step-function
plot_bas_hazard(result, xlab = 'Time [semesters]')

# Plot frailty standard deviation
plot_frailty_sd(result, ylim=c(0, 1), xlab = 'Time [intervals]', ylab = 'Standard deviation')

# Plot posterior frailty estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
plot_post_frailty_est(result, data_dropout,
                      xlim = c(1, 12), ylim = c(0, 3),
                      pch_type = pch_type, color_bg = color_bg,
                      xlab = 'Time [intervals]', ylab = 'Posterior estimates',
                      pos_legend = 'bottomright')

#-------------------------------------------------------------------------------
# Compute the frailty standard deviation in the reduced mode
reduced_frailty_sd <- frailty_sd.AdPaik(result, FALSE)

# Plot frailty standard deviation
plot_frailty_sd(result, frailty_sd = reduced_frailty_sd, flag_variance = TRUE,
                ylim=c(0, 1), main_title = 'Frailty variance')

# Plot posterior frailty estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
plot_post_frailty_est(result, data_dropout, 
                      flag_eps = TRUE, flag_alpha = FALSE,
                      ylim=c(0.3,1), xlim =c (1, 12),
                      pch_type = pch_type, color_bg = color_bg,
                      xlab = 'Time [intervals]', ylab = 'Posterior estimates',
                      pos_legend = 'bottomright')

plot_post_frailty_est(result, data_dropout,
                      flag_eps = FALSE, flag_alpha = TRUE,
                      xlim=c(0.95, 1.05), ylim = c(0.3,1.05), 
                      pch_type = pch_type, color_bg = color_bg,
                      xlab = 'Time [intervals]', ylab = 'Posterior estimates',
                      pos_legend = 'bottomright')

#-------------------------------------------------------------------------------
# One dimensional analysis of the log-likelihood function
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)

index_param_to_vary <- 14
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = TRUE, optimal_params = result$OptimalParameters,
                             categories_range_min, categories_range_max, n_iter = 1)
analysis_1D_opt$EstimatedParameter
analysis_1D_opt$OptimizedLoglikelihood

result$OptimalParameters[index_param_to_vary]


analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = TRUE, optimal_params = result$OptimalParameters,
                             categories_range_min, categories_range_max, n_iter = 1)


eps = eps_paik
categories_range_min <- c(-8, -1, eps, eps, eps)
categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)
index_param_to_vary <- 12
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = FALSE, optimal_params = NULL,
                             categories_range_min, categories_range_max, n_iter = 5, flag_plot = TRUE)
analysis_1D_opt

