# Check Existence of Provided Input Index

The method controls that the provided input index exists: it cannot be
greater than the maximum number of parameters of the current model.

## Usage

``` r
check.index(index, n_params)
```

## Arguments

- index:

  Index with respect to which the user wants to study the one
  dimensional behaviour of the log-likelihood function.

- n_params:

  Number of parameters of the model

## Value

An error if any condition is not satisfied.
