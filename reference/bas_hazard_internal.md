# Internal Function for the Baseline Hazard Step-Function

The method computes the baseline hazard step-function in each interval
of the time-domain, using the estimated parameters \\\phi_k, \forall k\\

## Usage

``` r
bas_hazard_internal(optimal_params, time_axis)
```

## Arguments

- optimal_params:

  Numerical vector of length equal to the number of model parameters,
  containing the optimal estimated parameters.

- time_axis:

  Numerical vector of temporal domain.

## Value

Numerical vector of length equal to the number of intervals of the
time-domain, with the value of the baseline hazard step-function.
