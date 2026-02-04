# Evaluation of Model Group log-Likelihood

Evaluation of model group log-likelihood at the provided parameter
vector and data. This function is internally called by 'll_AdPaik_eval'
to evaluate the log-likelihood function, considering all and only the
individuals belonging to a group.

## Usage

``` r
ll_AdPaik_centre_eval(params, dataset, dropout_matrix, e_matrix)
```

## Arguments

- params:

  Parameter vector.

- dataset:

  Matrix of dataset regressors, with a number of rows equal to the
  number of individuals in a cluster.

- dropout_matrix:

  Binary matrix indicating in which interval of the time domain and
  individual failed. For an individual, the sum of the row elements must
  be equal to 1 (if he/she failed) or 0 (if he/she does not failed). It
  has dimension equal to (n_individuals, n_intervals)

- e_matrix:

  Matrix of dimension (n_individual, n_intervals), where each element
  contains the evaluation of the temporal integral, performed through
  the function @time_int_eval.

## Value

Group log-likelihood evaluation
