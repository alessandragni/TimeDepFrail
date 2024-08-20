# MAC SETUP
setwd("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R")
load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/Functions.RData")
load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/data_dropout.RData")
save.image("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/Functions.RData")
#load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/dataless_time_varying_year2012.RData")

#-------------------------------------------------------------------------------
# WRONG INPUT VARIABLES ONLY FOR THE ADAPTED PAIK ET AL.'S MODEL
#-------------------------------------------------------------------------------
# Remove or add one category: ERROR
# Add one positive random element to the time domain vector: ERROR
# Negative element in time domain vector: ERROR
# Not-existent covariate in the dataset, but indicated in the formula: ERROR
# Not provided cluster variable: ERROR
# Not provided response: ERROR (but ~ must be present, otherwise it complains.)
# Remove one variable in the model call: ERROR (at least one variable is mssing, with no default)
# Null value in the dataset: ERROR
# Only one center in the centre variable: ERROR

eps_paik <- 1e-10
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)

time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)

formula <- time_to_event ~ Gender + CFUP + cluster(centre)

result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)
#-------------------------------------------------------------------------------
plot.bas_hazard(result, xlim=c(1,result$TimeDomain[result$NIntervals+1]))

plot.frailty_sd(result, ylim=c(0, 0.80))

pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))
plot.post_frailty_est(result, data_dropout,
                      ylim=c(0,2),
                      pch_type = pch_type, color_bg = color_bg)

reduced_frailty_sd <- frailty.sd.AdPaik(result, FALSE)

plot.frailty_sd(result, frailty_sd = reduced_frailty_sd, flag_variance = TRUE,
                ylim=c(0, 0.40), main_title = 'Frailty variance')
#-------------------------------------------------------------------------------
index_param_to_vary <- 2
analysis_1D_var <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = FALSE, optimal_params = NULL,
                             categories_range_min, categories_range_max)

index_param_to_vary <- 2
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = TRUE, optimal_params = result$OptimalParameters,
                             categories_range_min, categories_range_max)


