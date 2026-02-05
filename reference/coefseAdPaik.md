# Extracts the Standard Errors of the Coefficients for the 'Adapted Paik et Al.' Model

Extracts the standard errors for \\\boldsymbol{\beta}\\ obtained with
the time-dependent frailty model proposed in the 'Adapted Paik et al.'
framework.

## Usage

``` r
coefseAdPaik(object)
```

## Arguments

- object:

  An S3 object of class `AdPaik`, returned by the main model function
  (`AdPaikModel`). This object contains all the optimal parameter
  estimates.

## Value

A named list containing the categories of the standard errors for the
optimal parameters.

## Details

The `se.coef` function extracts the standard errors for the estimated
parameters from the `StandardErrorParameters` field in `object`.

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
coefseAdPaik(result)
#> GenderMale       CFUP 
#> 0.05241492 0.03447011 
# }
```
