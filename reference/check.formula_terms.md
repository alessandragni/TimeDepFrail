# Check Correctness of Formula Terms

The function controls that the terms composing the formula object
provided in input to the main model are correct. They must include:

- response variable on the left hand side

- covariates (numerical or categorical) on the right hand side

- cluster variable (categorical) on the right hand side and specified by
  the function 'cluster()'

Moreover, it controls that the covariates are contained in the dataset
provided.

## Usage

``` r
check.formula_terms(formula, data)
```

## Arguments

- formula:

  Formula object specifying the relationship between the time-to-event,
  the covariates and the cluster variables.

- data:

  Dataset in which these variables can be found.

## Value

An error if any condition is not satified.
