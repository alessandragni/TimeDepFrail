# Evaluation of Model log-Likelihood

Evaluation of the log-likelihood function at the provided parameter
vector and data.

## Usage

``` r
ll_AdPaik_eval(params, dataset, centre, time_axis, dropout_matrix, e_matrix)
```

## Arguments

- params:

  Parameter vector

- dataset:

  Matrix of dimension equal to (number of individuals in the study,
  number of regressors), where only the regressors indicated in the
  formula object are considered.

- centre:

  Vector of length equal to the number of individuals in the study,
  where each element corresponds to the individual cluster membership.

- time_axis:

  Temporal domain

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

Overall log-likelihood function value at the provided parameters and
data

## Details

The function divides the individuals according to their group/cluster
membership and then evaluates the group log-likelihood through another
implemented function, but using all and only the individuals belonging
to that group. Then the results are summed together to return the
overall log-likelihood value.
