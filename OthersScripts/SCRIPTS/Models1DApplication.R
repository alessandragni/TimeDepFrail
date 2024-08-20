# WINDOWS SETUP
setwd("C:/Users/admin/Documents/R/Thesis/TimeVaryingSharedFrailtyCoxModels-R")
load("C:/Users/admin/Documents/R/Thesis/TimeVaryingSharedFrailtyCoxModels-R/Data/dataless_time_varying_year2012.RData") 

#-------------------------------------------------------------------------------
# MODEL APPLICATION
dataset_2012 <- as.matrix(data_app_year[,1:2])
time_to_event_2012 <- data_app_year[,3]
centre <- faculty_codes
time_axis <- a_interval

# Dataframe
dataset_2012 <- data.frame(dataset_2012)
dataset_2012 <- cbind(dataset_2012, time_to_event_2012, centre)
colnames(dataset_2012) <- c("GenderF", "CFUP", "time_to_event", "centre")
#-------------------------------------------------------------------------------
## ADAPTED PAIK ET AL MODEL
# Model call
eps_paik <- 1e-10
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0, 1 - eps_paik, 1, 10)

formula <- time_to_event ~ GenderF + CFUP + cluster(centre)

result <- AdPaikModel(formula, dataset_2012, time_axis, 
                      categories_range_min, categories_range_max)
summary(result)

# Perform one-dimensional analysis with respect to indicated parameter
index_param_to_vary <- 14

# With optimal parameters
AdPaik_1D(formula, dataset_2012, time_axis, index_param_to_vary, 
          result$OptimalParameters, flag_optimal_params = TRUE,
          categories_range_min, categories_range_max,
          n_iter = 5, tol_optimize = 1e-6,
          flag_plot = TRUE)

# With no optimal parameters
AdPaik_1D(formula, dataset_2012, time_axis, index_param_to_vary,
          0, flag_optimal_params = FALSE,
          categories_range_min, categories_range_max,
          n_iter = 5, tol_optimize = 1e-6,
          flag_plot = TRUE)

#-------------------------------------------------------------------------------
# CENTRE-SPECIFIC FRAILTY MODEL WITH POWER PARAMETER
eps_pp <- 1e-10
categories_range_min <- c(-8, -2, eps_pp, eps_pp)
categories_range_max <- c(-eps_pp, 0, 10, 1 - eps_pp)

C_mult <- 1

formula <- time_to_event ~ GenderF + CFUP + cluster(centre)

result <- PowParModel(formula, dataset_2012, time_axis, 
                      categories_range_min, categories_range_max, C_mult)
summary(result)

# Perform one-dimensional analysis with respect to indicated parameter
index_param_to_vary <- 5

# With optimal parameters
PowPar_1D(formula, dataset_2012, time_axis, index_param_to_vary, 
          result$OptimalParameters, flag_optimal_params = TRUE,
          categories_range_min, categories_range_max, C_mult,
          n_iter = 5, tol_optimize = 1e-6,
          flag_plot = TRUE)

# With no optimal parameters
PowPar_1D(formula, dataset_2012, time_axis, index_param_to_vary,
          0, flag_optimal_params = FALSE,
          categories_range_min, categories_range_max, C_mult,
          n_iter = 5, tol_optimize = 1e-6,
          flag_plot = TRUE)

#-------------------------------------------------------------------------------
# STOCHASTIC TIME-DEPENDENT CENTRE-SPECIFIC FRAILTY MODELS
eps_log <- 1e-6
categories_range_min <- c(-8, -2, eps_log, eps_log, eps_log)
categories_range_max <- c(-eps_log, 2, 2, 2, pi)

C_mult <- 1

formula <- time_to_event ~ GenderF + CFUP + cluster(centre)

result <- StocTimeDepModel(formula, dataset_2012, time_axis, 
                           categories_range_min, categories_range_max, C_mult)
summary(result)

# Perform one-dimensional analysis with respect to indicated parameter
index_param_to_vary <- 5

# With optimal parameters
StocTimeDep_1D(formula, dataset_2012, time_axis, index_param_to_vary, 
               result$OptimalParameters, flag_optimal_params = TRUE,
               categories_range_min, categories_range_max, C_mult,
               n_iter = 5, tol_optimize = 1e-6,
               flag_plot = TRUE)

# With no optimal parameters
StocTimeDep_1D(formula, dataset_2012, time_axis, index_param_to_vary,
               0, flag_optimal_params = FALSE,
               categories_range_min, categories_range_max, C_mult,
               n_iter = 2, tol_optimize = 1e-6,
               flag_plot = TRUE)
