#-------------------------------------------------------------------------------
#' @title
#' Internal function to check the structure of the model output
#'
#' @description
#' This internal function checks if the input `result` object belongs to one of the 
#' expected S3 classes ('AdPaik', 'PowPar', or 'StocTimeDep') and calls the 
#' appropriate checking function based on the class.
#'
#' @param result S3 object. Expected to be of class 'AdPaik', 'PowPar', or 'StocTimeDep'.
#'
#' @details
#' This function is internal and dispatches the appropriate structure-checking function 
#' based on the class of the input. Currently, it supports 'AdPaik' class. The structure-checking 
#' functions for 'PowPar' and 'StocTimeDep' are placeholders and should be implemented if needed.
#'
#' @return Throws an error if the `result` object does not belong to one of the expected classes, 
#' or if the object’s structure is incorrect.
#'
#' @keywords internal

check.result <- function(result){
  # Check the class of result
  if(!inherits(result, 'AdPaik')) # & (!inherits(result, 'PowPar')) & (!inherits(result, 'StocTimeDep')))
    stop("Wrong S3 class object for input 'result' argument.")
  
  if(inherits(result, 'AdPaik'))
    check.result.AdPaik(result)
  #else if(inherits(result, 'PowPar'))
  #  check.result.PowPar(result)
  #else
  #  check.result.StocTimeDep(result)
}

#-------------------------------------------------------------------------------
#' @title
#' Check structure of the 'AdPaikModel' output
#'
#' @description
#' The function controls that the structure of the input variable is coherent with the one returned by the
#' 'AdPaikModel' execution.
#'
#' @param result S3 object of class 'AdPaik', composed of several elements. See details.
#'
#' @details The output of the model call 'AdPaikModel(...)' is a S3 object of class 'AdPaik', composed of 18 quantities:
#' - formula: formula object provided in input by the user and specifying the relationship between the time-to-event, the covariates of
#' the dataset (regressors) and the cluster variable.
#' - Regressors: categorical vector of length R, with the name of the regressors.
#' They could be different from the original covariates of the dataset in case of categorical covariates.
#' Indeed, each categorical covariate with n levels needs to be transformed into (n-1) dummy variables and, therefore, (n-1) regressors.
#' - NRegressors: number of regressors (R)
#' - ClusterVariable: name of the variable with respect to which the individuals can be grouped.
#' - NClusters: number of clusters/groups/centres
#' - NIntervals: number of intervals of the time-domain, also called with L. It corresponds to the length of the time-domain minus 1.
#' - NParameters: number of parameters of the model. It can be computed as: \eqn{n_p = 2L + R + 2}.
#' - ParametersCategories: Numerical vector of length 5, containing the numerosity of each parameter category.
#' - ParametersRangeMin: Numerical vector of length \eqn{n_p}, giving the minimum range of each parameter.
#' - ParametersRangeMax: Numerical vector of length \eqn{n_p}, giving the maximum range of each parameter.
#' - Loglikelihood: value of the maximized log-likelihood function, at the optimal estimated parameters.
#' - AIC: 'Akaike Information Criterion': it can be computed as \eqn{AIC = 2n_p - 2ll_{optimal}}.
#' It gives an idea of the loss of information related to the model fitting and output.
#' The smaller it is, the less loss of information and the better model accuracy.
#' - Status: Logical value. Does the model reach convergence? If so, the variable is TRUE, otherwise FALSE.
#' - NRun: Number of runs necessary to reach convergence. If the model does not reach convergence, such number is equal to the maximum number
#' of imposed runs.
#' - OptimalParameters: numerical vector of length \eqn{n_p}, containing the optimal estimated parameters or, in other words, the parameters
#' that maximizes the log-likelihood function.
#' - StandardErrorParameters: numerical vector of length \eqn{n_p}, corresponding to the standard error of each estimated parameters.
#' - ParametersCI: S3 object of class 'ParametersCI', composed of two numerical vector of length equal to \eqn{n_p}: the left and right confidence
#' interval of each estimated parameter.
#' - FrailtyStandardDeviation: numerical vector of length equal to L (i.e. number of intervals of the time-domain), reporting the standard deviation
#' of the frailty.
#' - PosteriorFrailtyEstimates: S3 object of class 'PFE.AdPaik'. See details.
#' - PosteriorFrailtyVariance: S3 object of class 'PFV.AdPaik'. See details.
#' - PosteriorFrailtyCI: S3 object of class 'PFCI.AdPaik'. See details.
#'
#' @details
#' The object of class 'PFE.AdPaik' contains the Posterior Frailty Estimates computed with the procedure indicated in the reference paper and
#' it is composed of three elements:
#' - 'alpha': posterior frailty estimates for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of groups/centres.
#' - 'eps': posterior frailty estimates for \eqn{\epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - 'Z': posterior frailty estimates for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'
#' @details
#' The object of class 'PFV.AdPaik' contains the Posterior Frailty Variances computed as indicated in the reference papaer and it
#' is  composed of three elements:
#' - 'alphaVar': posterior frailty variance for \eqn{\alpha_j, \forall j}. It is a vector of length equal to the number of groups/centres.
#' - 'epsVar': posterior frailty variance for \eqn{\epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - 'ZVar': posterior frailty variance for \eqn{Z_{jk} = \alpha_j + \epsilon_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'
#' @details
#' The object of class 'PFCI.AdPaik' contains the Posterior Frailty Confidence Interval and it is composed of two elements:
#' - left confidence interval for the estimated \eqn{\hat{Z}_{jk}, \forall j,k}. Matrix of dimension (N, L).
#' - right confidence interval for the estimated \eqn{\hat{Z}_{jk}, \forall j,k}. Matrix of dimension (N, L).
#'
#'
#' @return An error if any condition is not satisfied.

check.result.AdPaik <- function(result){
  # Save the names of the list elements
  names_list.AdPaik <- c("formula", "Regressors", "NRegressors", "ClusterVariable", "NClusters",
                         "TimeDomain", "NIntervals",
                         "NParameters", "ParametersCategories",
                         "ParametersRange",
                         "Loglikelihood", "AIC", "Status", "NRun",
                         "OptimalParameters", "StandardErrorParameters",
                         "ParametersCI", "BaselineHazard",
                         "FrailtyDispersion", "PosteriorFrailtyEstimates",
                         "PosteriorFrailtyVariance","PosteriorFrailtyCI")
  
  # Other than a class, it is a list
  if(! is.list(result))
    stop("Wrong structure for input 'result' argument.")
  
  names_list <- names_list.AdPaik
  for(i in 1:length(names_list)){
    if(names(result)[i] != names_list[i]){
      msg <- paste(names_list[i], "does not appear in the input 'result' argument. ")
      stop(msg)
    }
  }
  
  # Compute the number of parameters
  n_params <- result$NIntervals * 2 + result$NRegressors + 2
  
  # For each element of the list, control its structure
  if(!inherits(result$formula, "formula"))
    stop("'formula' is not a formula object.")
  if(!is.vector(result$Regressors))
    stop("'Regressors' is not a vector.")
  if(!is.numeric(result$NRegressors))
    stop("'NRegressors' is not a number.")
  if(!is.character(result$ClusterVariable))
    stop("'ClusterVariable' is not a string.")
  if(!is.numeric(result$NClusters))
    stop("'NCluster' is not a number.")
  
  if(! is.vector(result$TimeDomain))
    stop("'TimeDomain' is not a vector.")
  if(!is.numeric(result$NIntervals))
    stop("'NIntervals' is not a number.")
  if(length(result$TimeDomain) - 1 != result$NIntervals)
    stop("Different values for number of intervals in 'TimeDomain' and 'NIntervals'")
  
  if(! is.numeric(result$NParameters))
    stop("'NParameters' is not a number.")
  if(! is.vector(result$ParametersCategories))
    stop("'ParametersCategories' is not a vector.")
  if(length(result$ParametersCategories) != 5)
    stop("Wrong length of 'ParametersCategories' vector.")
  
  check.params_range(result$ParametersRange, n_params)
  
  if(! is.numeric(result$Loglikelihood))
    stop("'Loglikelihood' is not a value.")
  if(! is.numeric(result$AIC))
    stop("'AIC' is not a value.")
  if(! is.logical(result$Status))
    stop("'Status' is not a binary variable.")
  if(! is.numeric(result$NRun))
    stop("NRun' is not a number.")
  
  if(! is.vector(result$OptimalParameters))
    stop("'OptimalParameters' is not a vector.")
  if(! is.vector(result$StandardErrorParameters))
    stop("'StandardErrorParameters' is not a vector.")
  if(length(result$OptimalParameters) != n_params)
    stop("Wrong length of 'OptimalParameters' vector.")
  if(length(result$StandardErrorParameters) != n_params)
    stop("Wrong length of 'StandardErrorParameters' vector.")
  
  check.structure_paramsCI(result$ParametersCI)
  
  if(! is.vector(result$BaselineHazard))
    stop("'BaselineHazard' is not a vector.")
  if(length(result$BaselineHazard) != result$NIntervals)
    stop("Wrong length of 'BaselineHazrad' vector.")
  
  check.frailty_dispersion(result$FrailtyDispersion, result$NIntervals)
  check.pos_frailty_sd(result$FrailtyDispersion$FrailtyStandardDeviation, result$NIntervals)
  
  check.structure_post_frailty_est(result$PosteriorFrailtyEstimates, result$NIntervals, result$NClusters)
  check.structure_post_frailty_var(result$PosteriorFrailtyVariance, result$NIntervals, result$NClusters)
  check.structure_post_frailty_CI(result$PosteriorFrailtyCI, result$NIntervals, result$NClusters)
}

