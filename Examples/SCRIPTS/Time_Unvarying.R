
library(TimeDepFrail)
data("data_dropout")

data_dropout2 = data_dropout
data_dropout2$status = ifelse(data_dropout2$time_to_event < 6.1, 1, 0)


library(frailtypack)

# SHARED FRAILTY COX HAZARD MODEL
frailty_model <- frailtyPenal(Surv(time_to_event, status) ~ cluster(group) + Gender + CFUP, data = data_dropout2, 
                              cross.validation = FALSE,
                              n.knots = 15, kappa = 5, hazard="Splines")
summary(frailty_model)

# Regression coefficients and their p_values
coeff_reg <- frailty_model$coef
hazard_ratio_frailty <- exp(coeff_reg)
p_values_coeff_reg <- frailty_model$beta_p.value

# # Comparison with Cox-PH base
# summary(cox_phmodel)
# coeff_reg_cox_base <- cox_phmodel$coefficients

# Variance of the frailty parameter
theta <- frailty_model$theta
p_value_theta <- frailty_model$theta_p.value

# Marginal penalized log-likelihood in the semiparametric case
ploglik <- frailty_model$logLikPenal

# Number of observations and of events
n_events <- frailty_model$n.events
n_obs <- frailty_model$n
n_iter <- frailty_model$n.iter
#n_group <- frailty_model$groups

#-------------------------------------------------------------------------------
# Baseline Survival and Hazard function 
frailty_time <- frailty_model$x
frailty_hazard <- frailty_model$lam[,1,1]
frailty_survival <- frailty_model$surv[,1,1]

# Plot with default commands on the whole time domain
plot(frailty_model, type.plot = "Survival", conf.bands=FALSE,
     main = "Estimated baseline survival function", color="darkorange", median=TRUE, Xlab = "Time [semesters]", Ylab = "Survival probability")

# dev.new()
# plot(frailty_model, type.plot = "Hazard", conf.bands=FALSE,
#      main = "Estimated baseline Hazard function", color="black", median=TRUE, Xlab = "Time[semester]", Ylab = "Instantaneous risk of failure")

# Plot corrected baseline hazard function (area divided)
smoothingSpline = smooth.spline(frailty_time,frailty_hazard, spar=0.35)
area<-0 # normalize Spline
for(ii in 1:(length(smoothingSpline$y)-1)){
  area<-area+(smoothingSpline$x[ii+1]-smoothingSpline$x[ii])*smoothingSpline$y[ii+1]
}
smoothingSpline$y<-smoothingSpline$y/area
plot(smoothingSpline, type="l",
     main = "Estimated baseline hazard function", color="black", median=TRUE, xlab = "Time [semesters]", ylab = "Instantaneous risk of failure")

plot(smoothingSpline$x[17:100], smoothingSpline$y[17:100], type = 'l', col="black",
     main = "Estimated baseline hazard function", xlab = "Time [semesters]", ylab="Instantaneous risk of failure")
