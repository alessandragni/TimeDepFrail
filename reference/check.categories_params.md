# Check Correctness of Parameters Categories

The function controls that the provided parameters categories have a
length equal to the number of categories required by the model
parameters. For the current model, the number of categories is 5 because
there are five blocks of unkown parameters (\\\phi_k \forall k, \beta_r
\forall r, \mu_1, \nu, \gamma_k \forall k\\).

Moreover, it also controls that the minimum value of a parameter
category is actually less than or eqaul to the maximum value for the
same category.

## Usage

``` r
check.categories_params(
  n_categories,
  categories_range_min,
  categories_range_max
)
```

## Arguments

- n_categories:

  Number of categories expected by the model. For the current model,
  they are 5.

- categories_range_min:

  Numerical vector of length 5, containing the minimum ranges for the
  parameters beloning to those categories.

- categories_range_max:

  Numerical vector of length equal to 5, containing the maximum ranges
  for the parameters belonging to those categories.

## Value

An error if the any condition is not satisfied.
