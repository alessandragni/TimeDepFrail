# Extract AIC for `AdPaik` Objects

Computes the AIC for an `AdPaik` model.

## Usage

``` r
# S3 method for class 'AdPaik'
extractAIC(fit, scale = NULL, k = 2, ...)
```

## Arguments

- fit:

  An `AdPaik` model object.

- scale:

  Changing it is not supported for this model. It will be ignored.

- k:

  Penalty parameter (default is 2 for AIC).

- ...:

  Additional arguments (ignored).

## Value

A numeric vector with the number of parameters and AIC value.
