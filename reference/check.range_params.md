# Check Correctness of Input Parameters

The function controls that the input parameter vector have a length
equal to the theoretical one required by the model and that each
parameter properly belongs to its range.

## Usage

``` r
check.range_params(optimal_params, params_range_min, params_range_max)
```

## Arguments

- optimal_params:

  Numerical vector of length equal to the number of model parameters.
  For the 'Adapted Paik et al.'s Model' it can be computed as: \\n_p =
  2L + R + 2\\, where \\L\\ stands for the number of intervals of the
  time domain and \\R\\ the number of regressors of the dataset.

- params_range_min:

  Numerical vector of length equal to the number of model parameters
  (\\n_p\\) and containing the minimum range for each parameter.

- params_range_max:

  Numerical vector of length equal to the number of model parameters
  (\\n_p\\) and containing the maximum range for each parameter.

## Value

An error if any condition is not satisfied.
