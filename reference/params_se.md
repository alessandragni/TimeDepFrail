# Standard Error of the Parameters

Function for computing the standard error of each optimal parameter,
estimated through the constraint multi-dimensional optimization. The
procedure for the computation is based on the numerical approximation of
the second derivative of the log-likelihood function, by the 'centered
finite difference scheme' with an accuracy of the second order.

## Usage

``` r
params_se(
  optimal_params,
  params_range_min,
  params_range_max,
  dataset,
  centre,
  time_axis,
  dropout_matrix,
  e_matrix,
  h_dd
)
```

## Arguments

- optimal_params:

  Numerical vector of optimal parameters. Its length (i.e. number of
  parameters) is equal to \\n_p\\.

- params_range_min:

  Numerical vector of length equal to \\n_p\\, containing the minimum
  range of each parameter.

- params_range_max:

  Numerical vector of length equal to \\n_p\\, containing the maximum
  range of each parameter.

- dataset:

  Dataset containing the value of the regressors for all individuals in
  the study.

- centre:

  vector containing the group membership of each individual and that
  induces the clustering subdivision.

- time_axis:

  Temporal domain. Its number of intervals corresponds to the length of
  the time-domain minus 1

- dropout_matrix:

  Binary matrix of dimension (n_individuals, n_intervals). The sum of
  the elements of each row must be (1), if the associated individual
  failed in a precise interval, and (0) if the individual did not fail
  in the @time-axis. Therefore, if an individual failed in the
  time-domain, the interval in which he failed will have value (1) and
  the others (0).

- e_matrix:

  Matrix of dimension (n_individuals, n_intervals) where each element
  contains the resolution of the temporal integral for that individual
  in that interval, thorugh the 'e_time_fun' function.

- h_dd:

  Discretization step for the numerical approximation of the second
  derivative fo the loglikelihood function.

## Value

Vector of parameter standard error, of length equal to the number of
model parameters.

## Details

The standard error of each parameter is computed as the inverse of the
square root of the 'Information matrix', that in turn is computed as the
opposite of the 'Hessian matrix'. Only its diagonal is built and its
elements are separatey evaluated through a numerical approximation of
the second derivative of the log-likelihood function.

The function requires the optimal parameter vector and other
parameters-related variables, to check:

- the right numerosity of the parameter vector

- the correct range existence of each parameter (i.e. each parameter
  lies in its range).
