# Extracts the Confidence Intervals for the Coefficients for the 'Adapted Paik et Al.' Model

Extracts the confidence intervals for \\\boldsymbol{\beta}\\ obtained
with the time-dependent frailty model proposed in the 'Adapted Paik et
al.' framework.

## Usage

``` r
# S3 method for class 'AdPaik'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  An S3 object of class `AdPaik`, returned by the main model function
  (`AdPaikModel`). This object contains all the optimal parameter
  estimates.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. Defaults
  to NULL, and all parameters are considered. Changing it is not
  supported for this model. It will be ignored.

- level:

  The confidence level required. Defaults to 0.95.

- ...:

  Additional arguments to be passed to other methods.

## Value

A named list containing the categories of the standard errors for the
optimal parameters.

## Details

The `confint.AdPaik` function extracts the standard errors for the beta
coefficients from the `ParametersCI` field in `object`.

The function validates the structure of `object` and ensures
compatibility with the expected model output. It throws an error if the
object is malformed or inconsistent.

## Examples

``` r
# Example using the 'Academic Dropout' dataset
data(data_dropout)

# Define the formula and time axis for the model
formula <- time_to_event ~ Gender + CFUP + cluster(group)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)

# \donttest{
# Run the main model
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max, TRUE)

# Extract the coefficients
confint(result)
#>                 2.5 %     97.5 %
#> GenderMale -0.1027323  0.1027309
#> CFUP       -1.3422364 -1.2071159
# }
```
