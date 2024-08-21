
devtools::install_github("alessandragni/TimeDepFrail")

library(TimeDepFrail)
data("data_dropout")

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

# Plot baseline hazard step-function
plot_bas_hazard(result, xlim=c(1,result$TimeDomain[result$NIntervals+1]),
                xlab = 'Time [semesters]', ylab = 'Baseline hazard')

# Plot frailty standard deviation
plot_frailty_sd(result, ylim=c(0, 0.50), xlab = 'Time [intervals]', ylab = 'Standard deviation')

# Plot posterior frailty estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
plot_post_frailty_est(result, data_dropout,
                      ylim=c(0,3), xlim =c (1, 11),
                      pch_type = pch_type, color_bg = color_bg,
                      xlab = 'Time [intervals]', ylab = 'Posterior estimates',
                      pos_legend = 'bottomright')

#-------------------------------------------------------------------------------
# Compute the frailty standard deviation in the reduced mode
reduced_frailty_sd <- frailty_sd.AdPaik(result, FALSE)

# Plot frailty standard deviation
plot_frailty_sd(result, frailty_sd = reduced_frailty_sd, flag_variance = TRUE,
                ylim=c(0, 0.40), main_title = 'Frailty variance')

# Plot posterior frailty estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
plot_post_frailty_est(result, data_dropout, flag_esp = TRUE, flag_alpha = FALSE,
                      ylim=c(0,3), xlim =c (1, 11),
                      pch_type = pch_type, color_bg = color_bg,
                      xlab = 'Time [intervals]', ylab = 'Posterior estimates',
                      pos_legend = 'bottomright')

plot_post_frailty_est(result, data_dropout,
                      flag_eps = TRUE, flag_alpha = FALSE,
                      ylim=c(0.6,2.5),
                      pch_type = pch_type, color_bg = color_bg)

#-------------------------------------------------------------------------------
# One dimensional analysis of the log-likelihood function
index_param_to_vary <- 2
analysis_1D_var <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = FALSE, optimal_params = NULL,
                             categories_range_min, categories_range_max)

index_param_to_vary <- 2
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = TRUE, optimal_params = result$OptimalParameters,
                             categories_range_min, categories_range_max)
analysis_1D_opt$EstimatedParameter
result$OptimalParameters[index_param_to_vary]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# CENTRE-SPECIFIC FRAILTY MODEL WITH POWER PARAMETER
eps_pp <- 1e-10
categories_range_min <- c(-8, -2, eps_pp, eps_pp)
categories_range_max <- c(-eps_pp, 0.5, 10, 1 - eps_pp)

C_mult <- 1

time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)

formula <- time_to_event ~ Gender + CFUP + cluster(group)

result <- PowParModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max, C_mult)
summary(result)

# Plot baseline hazard step-function
plot_bas_hazard(result, xlim=c(1,result$TimeDomain[result$NIntervals+1]))

# Plot frailty standard deviation
plot_frailty_sd(result, ylim=c(0, 0.80))


#-------------------------------------------------------------------------------
# STOCHASTIC TIME-DEPENDENT CENTRE-SPECIFIC FRAILTY MODELS
eps_log <- 1e-6
range_min_lf <- c(rep(-8, n_interval), rep(-2, n_regressor), eps_log, eps_log, eps_log)
range_max_lf <- c(rep(-eps_log, n_interval), rep(0.5, n_regressor), 2, 2, pi)

C_mult <- 1

time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)

formula <- time_to_event ~ Gender + CFUP + cluster(group)

result <- StocTimeDepModel(formula, data_dropout, time_axis,
                           categories_range_min, categories_range_max, C_mult)
summary(result)

#-------------------------------------------------------------------------------
# Model call with global method and specifying the name of desired model
result.AdPaik <- TimeDepFrailty(formula, data_dropout, time_axis,
                         categories_range_min, categories_range_max, model_type = 'AdPaikModel')

result.PowPar <- TimeDepFrailty(formula, data_dropout, time_axis,
                                categories_range_min, categories_range_max, C_mult = 1,
                                model_type = 'PowParModel')

#-------------------------------------------------------------------------------






