devtools::install_github("alessandragni/TimeDepFrail")

library(TimeDepFrail)

help(package = "TimeDepFrail")

#### Section 4: Syntax and Implementation details #####

##### Section 4.1 #####
# ?AdPaikModel
# ?summary.AdPaik
# ?coef.AdPaik
# ?coefse
# ?confint.AdPaik
# ?plot_bas_hazard

##### Section 4.2 #####
# ?frailty_sd
# ?plot_frailty_sd

##### Section 4.3 #####
# ?plot_post_frailty_est

##### Section 4.4 #####
# ?survival
# ?plot_survival

##### Section 4.5 #####
# ?AdPaik_1D


#### Section 5: Worked example #####

##### Section 5.1 #####
data(data_dropout)
head(data_dropout)

formula <- time_to_event ~ Gender + CFUP + cluster(group)

time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)

eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)

###### Section 5.2 ######
# Main model and summary
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)
summary(result)
print(summary(result))

# Baseline Hazard function
result$BaselineHazard

pdf("Examples/Plots/BaselineHazard.pdf", width=8, height=5)
plot_bas_hazard(result)
dev.off()

# Estimated regressors and frailty parameters
coef(result)
coefse(result)
confint(result)


##### Section 5.3 #####

# Frailty standard deviation
result$FrailtyDispersion$FrailtyVariance
result$FrailtyDispersion$FrailtyStandardDeviation

pdf("Examples/Plots/FrailtySD.pdf", width=8, height=5)
plot_frailty_sd(result)
dev.off()

plot_frailty_sd(result, flag_variance = TRUE)


red_frailty_sd <- frailty_sd(result, flag_fullsd = FALSE)
plot_frailty_sd(result, frailty_sd = red_frailty_sd, flag_variance = FALSE)


# Posterior Frailty Estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))

pdf("Examples/Plots/PostFrailtyEst.pdf", width=8, height=5)
plot_post_frailty_est(result, data_dropout,
                      pch_type = pch_type, color_bg = color_bg)
dev.off()

mean(result$PosteriorFrailtyEstimates$Z)

pdf("Examples/Plots/post_estimates_eps.pdf", width=5, height=5)
plot_post_frailty_est(result, data_dropout,
                      flag_eps = TRUE, flag_alpha = FALSE,
                      pch_type = pch_type, color_bg = color_bg)
dev.off()

pdf("Examples/Plots/post_estimates_alpha.pdf", width=5, height=5)
plot_post_frailty_est(result, data_dropout,
                      flag_eps = FALSE, flag_alpha = TRUE,
                      pch_type = pch_type, color_bg = color_bg)
dev.off()


##### Section 5.4 ##### 

survival_df = survival(result, data_dropout)

pdf("Examples/Plots/Survival.pdf", width=8, height=5)
plot_survival(result, survival_df)
dev.off()


##### Section 5.5 ##### 

# Identify a parameter existence range
set.seed(1)
index_param_to_vary <- 1
analysis_1D_opt <- AdPaik_1D(formula, data_dropout,
                             time_axis, index_param_to_vary, 
                             flag_optimal_params = FALSE, 
                             optimal_params = NULL,
                             flag_plot = TRUE,
                             categories_range_min, categories_range_max, 
                             n_iter = 5)

analysis_1D_opt



categories_range_min <- c(-8, -1, eps, eps, eps)
categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)

set.seed(1)
index_param_to_vary <- 12
analysis_1D_opt <- AdPaik_1D(formula, data_dropout, 
                             time_axis, index_param_to_vary, 
                             flag_optimal_params = FALSE, 
                             optimal_params = NULL,
                             flag_plot = TRUE,
                             categories_range_min, categories_range_max, 
                             n_iter = 5)
analysis_1D_opt



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

frailty_time <- frailty_model$x
frailty_hazard <- frailty_model$lam[,1,1]

smoothingSpline = smooth.spline(frailty_time, frailty_hazard, spar=0.35)
area = 0 
for(ii in 1:(length(smoothingSpline$y)-1)){
  area<-area+(smoothingSpline$x[ii+1]-smoothingSpline$x[ii])*smoothingSpline$y[ii+1]
}
smoothingSpline$y <- smoothingSpline$y/area

pdf("Examples/Plots/Baseline.pdf", width=10, height=5)
plot(smoothingSpline$x[17:98], smoothingSpline$y[17:98], type = 'l', col="black",
     main = "Estimated baseline hazard function", xlab = "Time [semesters]", ylab="Instantaneous risk of failure")
dev.off()



