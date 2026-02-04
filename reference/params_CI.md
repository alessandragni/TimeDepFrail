# Confidence Interval for the Optimal Estimated Parameters

The function provides the confidence interval for each estimated
parameter, using the standard error computed through another method and
provided as second argument to the current function.

## Usage

``` r
params_CI(optimal_params, se_params, level)
```

## Arguments

- optimal_params:

  Numerical vector of optimal estimated parameters. Its length is equal
  to the number of model parameters.

- se_params:

  Numerical vector containing the standard error associated to each
  estimated parameter.

- level:

  A numeric value representing the confidence level.

## Value

A S3 object of class 'ParametersCI', composed of two numerical vector of
length equal to the number of model parameters:

- ParamsCI_left: left confidence interval for each parameter

- ParamsCI_right: right confidence interval for each parameter
