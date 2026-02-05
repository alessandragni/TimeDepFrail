# Adapted Paik et Al.'s Model: Time-Dependent Shared Frailty Cox Model

Function for applying the 'Adapted Paik et al.'s Model', a Cox Model
with time-dependent frailty, shared by individuals belonging to the same
group/cluster.

To generate time-dependence, the temporal domain is divided into a
certain number of intervals. The model log-likelihood function depends
on a certain number of parameters and is maximized with respect to all
of them, using a reinterpretation of the 'Powell's method in
multidimension', that is a multi-dimensional optimization method based
on repeated one-dimensional optimization of the log-likelihood function
(with respect to one parameter at the time). In this context, the
one-dimensional optimization is performed through the 'optimize' R
function. For more information about the unknown model parameters, their
type and numerosity refer to Details.

Several quantities are estimated at the end of the optimization phase,
such as optimal parameters, baseline hazard, frailty dispersion
(standard deviation and variance), posterior frailty estimates, with
their variance and confidence interval, conditional survival function,
Akaike Information Criterion (AIC), ...

## Usage

``` r
AdPaikModel(
  formula,
  data,
  time_axis,
  categories_range_min,
  categories_range_max,
  n_extrarun = 60,
  tol_ll = 1e-06,
  tol_optimize = 1e-06,
  h_dd = 0.001,
  verbose = FALSE,
  print_previous_ll_values = c(TRUE, 3)
)
```

## Arguments

- formula:

  Formula object having on the left hand side the time-to-event
  variable, that is the time-instant in which the individual failed. On
  the right hand side, it has the regressors and the cluster variable.

- data:

  Dataset in which all variables of the formula object must be found and
  contained. This dataset can also contain other variables, but they
  will not be considered. It can be either a dataframe or a matrix, but
  in both cases the name of each column must correspond to the name of
  the formula variables. It is not necessary to attach it (in case of a
  data.frame)

- time_axis:

  Temporal domain

- categories_range_min:

  Vector containing the minimum value assumable by each parameter
  category.

- categories_range_max:

  Vector containing the maximum value assumable by each parameter
  category.

- n_extrarun:

  Total maximum number of runs (iterations) are obtained summing
  n_extrarun to the number of parameters.

- tol_ll:

  Tolerance on the log-likelihood value.

- tol_optimize:

  Internal tolerance for the one-dimensional optimization through
  'optimize' R function.

- h_dd:

  Discretization step used for the numerical approximation of the second
  derivative of the log-likelihood function.

- verbose:

  Logical. If `TRUE`, detailed progress messages will be printed to the
  console. Defaults to `FALSE`.

- print_previous_ll_values:

  If we want to print the previous values of the log-likelihood
  function. This can be useful for controlling that the optimization
  procedure is proceeding in a monotone way and it does not oscillate.
  This argument is composed of two elements: TRUE/FALSE if we want or
  not to print the previous values and how many values we want to print
  on the console. Default is (TRUE, 3), so that only the previous 3
  values of the log-likelihood are printed.

## Value

S3 object of class 'AdPaik', composed of several elements. See Details.

## Details

Two observation needs to made about the time-domain:

- The time domain may coincide with the follow-up or it may be contained
  in it. Indeed, the left boundary can be greater than the beginning of
  the follow-up and, for instance, it can coincide with the
  time-instants in which the events begin to happen; conversely, the
  right boundary of the two must be the same.

- The partition of the time domain into intervals can be made according
  to two selected criteria: (1) using an already existent partition of
  the follow-up (2) using the shape of the baseline hazard function for
  a time independent model as reference: divide the time-domain
  according to regions in which it has a peak or a plateau.

The parameters with respect to which the log-likelihood function must be
optimized are:

- baseline log-hazard (number of parameters = number of intervals of the
  time-domain)

- data regressors

- \\\mu_1\\, \\\nu\\: parameters of the gamma distribution of
  \\\alpha_j\\ (time-independent/constant) (2 parameters)

- \\\gamma_k\\: parameters of the gamma distribution of
  \\\epsilon\_{jk}\\ (time-dependent) (number of parameters = number of
  intervals) Another model parameter is \\\mu_2\\ and it is get imposing
  the constraint that \\\mu_1 + \mu_2 = 1\\. As it can be notice, some
  parameters can be grouped into the same category (regressors, baseline
  log-hazard and so on) and we can easily constraint them assigning each
  category both a minimum and maximum range. The vector is structured as
  follows: (baseline log-hazard, regressors, \\\mu_1\\, \\\nu\\,
  \\\gamma_k\\) with dimension (n_intervals, n_regressors, 1, 1,
  n_intervals).

The output of the model call 'AdPaikModel(...)' is a S3 object of class
'AdPaik', composed of the following quantities:

- formula: formula object provided in input by the user and specifying
  the relationship between the time-to-event, the covariates of the
  dataset (regressors) and the cluster variable.

- dataset: matrix of the dataset containing the regressors and the dummy
  variables of the categorical covariates.

- Regressors: categorical vector of length R, with the name of the
  regressors. They could be different from the original covariates of
  the dataset in case of categorical covariates. Indeed, each
  categorical covariate with n levels needs to be transformed into (n-1)
  dummy variables and, therefore, (n-1) new regressors.

- NRegressors: number of regressors (R)

- ClusterVariable: name of the variable with respect to which the
  individuals can be grouped.

- NClusters: number of clusters/centres (also indicated with N).

- ClusterCodes: vector of length equal to the number of clusters,
  containing the codes of the clusters.

- TimeDomain: vector of the time-domain partition.

- NIntervals: number of intervals of the time-domain, also called with
  L.

- NObservations: number of observations of the dataset.

- NParameters: number of parameters of the model. It can be computed as:
  \\n_p = 2L + R + 2\\.

- ParametersCategories: Numerical vector of length 5, containing the
  numerosity of each parameter category.

- ParametersRange: S3 object of class 'ParametersRange', containing
  ParametersRangeMin and ParametersRangeMax, two numerical vectors of
  length \\n_p\\, giving the minimum and the maximum range of each
  parameter, respectively.

- Loglikelihood: value of the maximized log-likelihood function, at the
  optimal estimated parameters.

- AIC: 'Akaike Information Criterion': it can be computed as \\AIC =
  2n_p - 2ll\_{optimal}\\. It quantifies the loss of information related
  to the model fitting and output. The smaller, the less the loss of
  information and the better the model accuracy.

- Status: Logical value. TRUE if the model reaches convergence, FALSE
  otherwise.

- NRun: Number of runs necessary to reach convergence. If the model does
  not reach convergence, such number is equal to the maximum number of
  imposed runs.

- OptimalParameters: numerical vector of length \\n_p\\, containing the
  optimal estimated parameters (the parameters that maximize the
  log-likelihood function).

- StandardErrorParameters: numerical vector of length \\n_p\\,
  corresponding to the standard error of each estimated parameters.

- ParametersCI: S3 object of class 'ParametersCI', composed of two
  numerical vector of length equal to \\n_p\\: the left and right 95\\
  interval of each estimated parameter of given level.

- BaselineHazard: numerical vector of length equal to L, containing the
  baseline hazard step-function.

- FrailtyDispersion: S3 object of class 'FrailtyDispersion', containing
  two numerical vectors of length equal to L with the standard deviation
  and the variance of the frailty. numerical vector of length equal to L
  (i.e. number of intervals of the time-domain), reporting the standard
  deviation of the frailty.

- PosteriorFrailtyEstimates: S3 object of class 'PFE.AdPaik'. The object
  of class 'PFE.AdPaik' contains the Posterior Frailty Estimates
  computed with the procedure indicated in the reference paper and it is
  composed of three elements:

  - 'alpha': posterior frailty estimates for \\\alpha_j, \forall j\\. It
    is a vector of length equal to the number of centres.

  - 'eps': posterior frailty estimates for \\\epsilon\_{jk}, \forall
    j,k\\. Matrix of dimension (N, L).

  - 'Z': posterior frailty estimates for \\Z\_{jk} = \alpha_j +
    \epsilon\_{jk}, \forall j,k\\. Matrix of dimension (N, L).

- PosteriorFrailtyVariance: S3 object of class 'PFV.AdPaik'. The object
  of class 'PFV.AdPaik' contains the Posterior Frailty Variances
  computed as indicated in the reference papaer and it is composed of
  three elements:

  - 'alphaVar': posterior frailty variance for \\\alpha_j, \forall j\\.
    It is a vector of length equal to the number of centres.

  - 'epsVar': posterior frailty variance for \\\epsilon\_{jk}, \forall
    j,k\\. Matrix of dimension (N, L).

  - 'ZVar': posterior frailty variance for \\Z\_{jk} = \alpha_j +
    \epsilon\_{jk}, \forall j,k\\. Matrix of dimension (N, L).

## Examples

``` r
# Consider the 'Academic Dropout dataset'
data(data_dropout)

# Define the variables needed for the model execution
formula <- time_to_event ~ Gender + CFUP + cluster(group)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)

# \donttest{
# Call the main model
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)
# }
```
