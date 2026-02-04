# Check Structure of the 'AdPaikModel' Output

The function controls that the structure of the input variable is
coherent with the one returned by the 'AdPaikModel' execution.

## Usage

``` r
check.result(result)
```

## Arguments

- result:

  S3 object of class 'AdPaik', composed of several elements. See
  details.

## Value

An error if any condition is not satisfied.

## Details

The output of the model call 'AdPaikModel(...)' is a S3 object of class
'AdPaik', composed of:

- formula: formula object provided in input by the user and specifying
  the relationship between the time-to-event, the covariates of the
  dataset (regressors) and the cluster variable.

- Regressors: categorical vector of length R, with the name of the
  regressors. They could be different from the original covariates of
  the dataset in case of categorical covariates. Indeed, each
  categorical covariate with n levels needs to be transformed into (n-1)
  dummy variables and, therefore, (n-1) regressors.

- NRegressors: number of regressors (R)

- ClusterVariable: name of the variable with respect to which the
  individuals can be grouped.

- NClusters: number of clusters/groups/centres

- ClusterCodes

- TimeDomain

- NIntervals: number of intervals of the time-domain, also called
  with L. It corresponds to the length of the time-domain minus 1.

- NParameters: number of parameters of the model. It can be computed as:
  \\n_p = 2L + R + 2\\.

- ParametersCategories: Numerical vector of length 5, containing the
  numerosity of each parameter category.

- ParametersRangeMin: Numerical vector of length \\n_p\\, giving the
  minimum range of each parameter.

- ParametersRangeMax: Numerical vector of length \\n_p\\, giving the
  maximum range of each parameter.

- Loglikelihood: value of the maximized log-likelihood function, at the
  optimal estimated parameters.

- AIC: 'Akaike Information Criterion': it can be computed as \\AIC =
  2n_p - 2ll\_{optimal}\\. It gives an idea of the loss of information
  related to the model fitting and output. The smaller it is, the less
  loss of information and the better model accuracy.

- Status: Logical value. Does the model reach convergence? If so, the
  variable is TRUE, otherwise FALSE.

- NRun: Number of runs necessary to reach convergence. If the model does
  not reach convergence, such number is equal to the maximum number of
  imposed runs.

- OptimalParameters: numerical vector of length \\n_p\\, containing the
  optimal estimated parameters or, in other words, the parameters that
  maximizes the log-likelihood function.

- StandardErrorParameters: numerical vector of length \\n_p\\,
  corresponding to the standard error of each estimated parameters.

- ParametersCI: S3 object of class 'ParametersCI', composed of two
  numerical vector of length equal to \\n_p\\: the left and right
  confidence interval of each estimated parameter.

- FrailtyStandardDeviation: numerical vector of length equal to L (i.e.
  number of intervals of the time-domain), reporting the standard
  deviation of the frailty.

- PosteriorFrailtyEstimates: S3 object of class 'PFE.AdPaik'. See
  details.

- PosteriorFrailtyVariance: S3 object of class 'PFV.AdPaik'. See
  details.

The object of class 'PFE.AdPaik' contains the Posterior Frailty
Estimates computed with the procedure indicated in the reference paper
and it is composed of three elements:

- 'alpha': posterior frailty estimates for \\\alpha_j, \forall j\\. It
  is a vector of length equal to the number of groups/centres.

- 'eps': posterior frailty estimates for \\\epsilon\_{jk}, \forall
  j,k\\. Matrix of dimension (N, L).

- 'Z': posterior frailty estimates for \\Z\_{jk} = \alpha_j +
  \epsilon\_{jk}, \forall j,k\\. Matrix of dimension (N, L).

The object of class 'PFV.AdPaik' contains the Posterior Frailty
Variances computed as indicated in the reference papaer and it is
composed of three elements:

- 'alphaVar': posterior frailty variance for \\\alpha_j, \forall j\\. It
  is a vector of length equal to the number of groups/centres.

- 'epsVar': posterior frailty variance for \\\epsilon\_{jk}, \forall
  j,k\\. Matrix of dimension (N, L).

- 'ZVar': posterior frailty variance for \\Z\_{jk} = \alpha_j +
  \epsilon\_{jk}, \forall j,k\\. Matrix of dimension (N, L).
