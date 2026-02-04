# Internal Function for Frailty Standard Deviation for the 'Adapted Paik et Al.'s Model'

The function computes both the standard deviation and the variance of
the time-dependent frailty of the 'Adapted Paik et al.'s Model'.

## Usage

``` r
frailty_Sd_internal(
  optimal_params,
  time_axis,
  n_regressors,
  categories_range_min,
  categories_range_max,
  flag_full
)
```

## Arguments

- optimal_params:

  Optimal parameter vector, estimated through multi-dimensional
  optimization of the log-likelihood function.

- time_axis:

  Partition of the temporal domain.

- n_regressors:

  Number of regressors of the dataset. This value must be provided to be
  able to correctly compute the number of parameters of the model and,
  therefore, to check the dimension of the parameter vector.

- categories_range_min:

  Vector of minimum value (range) assumed by the parameters category.

- categories_range_max:

  Vector of maximum value (range) assumed by the parameters category.

- flag_full:

  Do we want to compute the full frailty standard deviation (second
  case)? If so, the flag must be TRUE, otherwise (first case), FALSE.

## Value

S3 class object 'FrailtyDispersion' containing both two numerical
vectors of length equal to the number of intervals of the time-domain:

- FrailtyVariance

- FrailtyStandardDeviation
