# WINDOWS SETUP
# Set working directory, save the functions and load the dataset
setwd("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R")
load("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R/Data/Functions.RData")
load("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R/Data/data_dropout.RData")
#save.image("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R/Data/Functions.RData")
#load("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R/Data/dataless_time_varying_year2012.RData")

# MAC SETUP
setwd("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R")
load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/Functions.RData")
load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/data_dropout.RData")
#save.image("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/Functions.RData")
#load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/dataless_time_varying_year2012.RData")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# MODEL APPLICATION
#-------------------------------------------------------------------------------
# # Create dataset for model application: already saved in the workspace in the correct way
# data_dropout <- as.matrix(data_app_year[,1:2])
# time_to_event <- data_app_year[,3]
# centre <- faculty_codes
# time_axis <- a_interval
#
# # Dataframe
# data_dropout <- data.frame(data_dropout)
# data_dropout <- cbind(data_dropout, time_to_event, centre)
# colnames(data_dropout) <- c("Gender", "CFUP", "time_to_event", "group")
#-------------------------------------------------------------------------------
# RECALL TO CHANGE THE NUMBER OF RUNS IN THE FIRST TWO MODELS (internal variable)
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
plot.bas_hazard(result, xlim=c(1,result$TimeDomain[result$NIntervals+1]))

# Plot frailty standard deviation
plot.frailty_sd(result, ylim=c(0, 0.80))

# Plot posterior frailty estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
plot.post_frailty_est(result, data_dropout,
                      ylim=c(0,3),
                      pch_type = pch_type, color_bg = color_bg)

#-------------------------------------------------------------------------------
# Compute the frailty standard deviation in the reduced mode
reduced_frailty_sd <- frailty.sd.AdPaik(result, FALSE)

# Plot frailty standard deviation
plot.frailty_sd(result, frailty_sd = reduced_frailty_sd, flag_variance = TRUE,
                ylim=c(0, 0.40), main_title = 'Frailty variance')
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
plot.bas_hazard(result, xlim=c(1,result$TimeDomain[result$NIntervals+1]))

# Plot frailty standard deviation
plot.frailty_sd(result, ylim=c(0, 0.80))


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
save("C:/Users/admin/Documents/R/Thesis/TimeVaryingSharedFrailtyCoxModels-R/Data/Functions.RData")





