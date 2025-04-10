devtools::install_github("alessandragni/TimeDepFrail")

library(TimeDepFrail)

help(package = "TimeDepFrail")

# knitr::spin("ReplicationCode.R")

#### Section 5: Syntax and Implementation details #####

##### Section 5.1 #####
# ?AdPaikModel

# ?summary.AdPaik
# ?print.AdPaik
# ?plot.AdPaik
# ?logLik.AdPaik
# ?extractAIC.AdPaik
# ?nobs.AdPaik

# ?coef.AdPaik
# ?coefseAdPaik
# ?confint.AdPaik

# ?bas_hazard
# ?plot_bas_hazard

##### Section 5.2 #####
# ?frailty_sd
# ?plot_frailty_sd

# ?post_frailty_est
# ?plot_post_frailty_est
# ?post_frailty_var
# ?plot_post_frailty_var
# ?post_frailty_confint

##### Section 5.3 #####
# ?survivalAdPaik
# ?plot_survivalAdPaik

##### Section 5.4 #####
# ?AdPaik_1D

#### Section 6: Worked example #####

##### Section 6.1 #####
data(data_dropout)
head(data_dropout)

formula <- time_to_event ~ Gender + CFUP + cluster(group)

time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)

eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)

##### Section 6.2 #####
# Main model and summary
set.seed(1)
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)
summary(result)
print(result)
plot(result)
logLik(result)
extractAIC(result)
nobs(result)
formula(result)

# Estimated regressors
coef(result)
coefseAdPaik(result)
confint(result)

# Baseline Hazard function
bas_hazard(result)

pdf("Examples/Plots/BaselineHazard.pdf", width=8, height=5)
plot_bas_hazard(result)
dev.off()


##### Section 6.3 #####

# Frailty standard deviation and variance
frailty_sd(result)
frailty_sd(result, flag_variance = TRUE) 

pdf("Examples/Plots/FrailtySD.pdf", width=8, height=5)
plot_frailty_sd(result)
dev.off()

plot_frailty_sd(result, flag_variance = TRUE)

# Posterior Frailty Estimates

post_frailty_est(result, flag_eps = FALSE, flag_alpha = TRUE)
post_frailty_est(result, flag_eps = TRUE, flag_alpha = FALSE)
post_frailty_est(result, flag_eps = FALSE, flag_alpha = FALSE)

pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))

pdf("Examples/Plots/PostFrailtyEst.pdf", width=8, height=5)
plot_post_frailty_est(result, pch_type = pch_type, color_bg = color_bg)
dev.off()

mean(post_frailty_est(result))

pdf("Examples/Plots/post_estimates_eps.pdf", width=5, height=5)
plot_post_frailty_est(result, flag_eps = TRUE, flag_alpha = FALSE,
                      pch_type = pch_type, color_bg = color_bg)
dev.off()

pdf("Examples/Plots/post_estimates_alpha.pdf", width=5, height=5)
plot_post_frailty_est(result, flag_eps = FALSE, flag_alpha = TRUE,
                      pch_type = pch_type, color_bg = color_bg)
dev.off()

post_frailty_var(result)
plot_post_frailty_var(result)
post_frailty_confint(result)


##### Section 6.4 ##### 

survival_df = survivalAdPaik(result)

pdf("Examples/Plots/Survival.pdf", width=8, height=5)
plot_survivalAdPaik(result)
dev.off()


##### Section 6.5 ##### 

# Identify a parameter existence range
set.seed(123)
index_param_to_vary <- 1
analysis_1D_opt <- AdPaik_1D(formula, data_dropout,
                             time_axis, index_param_to_vary, 
                             flag_optimal_params = FALSE, 
                             optimal_params = NULL,
                             flag_plot = TRUE,
                             categories_range_min, categories_range_max, 
                             n_iter = 5)
analysis_1D_opt

# save the plots
pdf("Examples/Plots/ll_par1.pdf", width=8, height=5)
set.seed(123)
analysis_1D_opt <- AdPaik_1D(formula, data_dropout,
                             time_axis, index_param_to_vary, 
                             flag_optimal_params = FALSE, 
                             optimal_params = NULL,
                             flag_plot = TRUE,
                             categories_range_min, categories_range_max, 
                             n_iter = 5)
dev.off()


categories_range_min <- c(-8, -1, eps, eps, eps)
categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)

set.seed(123)
index_param_to_vary <- 12
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, 
                             time_axis, index_param_to_vary, 
                             flag_optimal_params = FALSE, 
                             optimal_params = NULL,
                             flag_plot = TRUE,
                             categories_range_min, categories_range_max, 
                             n_iter = 5)
analysis_1D_opt
# save the plots
pdf("Examples/Plots/ll_12par.pdf", width=8, height=5)
set.seed(123)
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, 
                             time_axis, index_param_to_vary, 
                             flag_optimal_params = FALSE, 
                             optimal_params = NULL,
                             flag_plot = TRUE,
                             categories_range_min, categories_range_max, 
                             n_iter = 5)
dev.off()


# Study the log-likelihood behaviour
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0.4, 1 - eps, 1, 10)
set.seed(1)
index_param_to_vary <- 14
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = TRUE, 
                             flag_plot = TRUE, optimal_params = result$OptimalParameters,
                             categories_range_min, categories_range_max, n_iter = 1)
analysis_1D_opt

# save the plots
pdf("Examples/Plots/Loglik14_1.pdf", width=8, height=5)
set.seed(1)
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, time_axis,
                             index_param_to_vary, flag_optimal_params = TRUE, 
                             flag_plot = TRUE, optimal_params = result$OptimalParameters,
                             categories_range_min, categories_range_max, n_iter = 1)
dev.off()



#### APPENDIX D ####

data("data_dropout")

data_dropout$status = ifelse(data_dropout$time_to_event < 6.1, 1, 0)
data_dropout$time_to_event = as.numeric(data_dropout$time_to_event)

library(frailtypack)

frailty_model <- frailtyPenal(Surv(time_to_event, status) ~ cluster(group) + Gender + CFUP, data = data_dropout, 
                              cross.validation = FALSE,
                              n.knots = 20, kappa = 1, hazard="Splines") 
# summary(frailty_model)

frailty_model$coef
frailty_model$beta_p.value


# Baseline Hazard function 

plot(frailty_model, type.plot = "Hazard", 
     xlim = c(1, 6.1), ylim = c(0, 0.045),
     main = "Estimated baseline hazard function", 
     col = 'black', conf.bands = FALSE)

pdf("Examples/Plots/Baseline.pdf", width=10, height=5)
plot(frailty_model, type.plot = "Hazard", 
     xlim = c(1, 6.1), ylim = c(0, 0.045),
     main = "Estimated baseline hazard function", 
     col = 'black', conf.bands = FALSE)
dev.off()

sessionInfo() 

