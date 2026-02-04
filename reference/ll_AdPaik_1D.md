# One-Dimensional Log-Likelihood Function to be Optimized

Model log-likelihood function to be optimized only with respect to a
parameter. To correctly identify this parameter inside the model and
inside the vector of all parameter, it is necessary to provide also the
position (index) of this parameter in the vector.

This function is internally used by the main function @AdPaikModel to
perform, as said, the one-dimensional optimization through 'optimize'.
It cannot be used to evaluate the log-likelihood function at a vector of
parameter and at the provided data. For this purpose, we have to use
another implemented function, called @ll_AdPaik_eval.

## Usage

``` r
ll_AdPaik_1D(
  x,
  index,
  params,
  dataset,
  centre,
  time_axis,
  dropout_matrix,
  e_matrix
)
```

## Arguments

- x:

  Value of the parameter, with respect to which the log-likelihood
  function has to be optimized.

- index:

  Index of the parameter inside the parameter vector. For instance, if
  we need to optimize the log-likelihood function with respect to the
  first regressor, then @x will be generic but @index will be equal to
  (n_intervals + 1) because in the parameter vector the first regressor
  appears after the baseline log-hazard group (n_intervals elements).

- params:

  Parameter vector.

- dataset:

  Matrix containing only the formula regressors, that is the regressors
  appearing in the formula object provided by the user and eventually
  modified if they are categorical (nd therefore transformed into dummy
  variables).

- centre:

  Individual membership to the clusters.

- time_axis:

  Temporal domain.

- dropout_matrix:

  Binary matrix indicating in which interval of the time domain an
  individual failed. For an individual, the sum of the row elements must
  be equal to 1 (if he/she failed) or 0 (if he/she does not failed). It
  has dimension equal to (n_individuals, n_intervals).

- e_matrix:

  Matrix of dimension (n_individual, n_intervals), where each element
  contains the evaluation of the temporal integral, performed through
  the function @param time_int_eval.

## Value

Overall log-likelihood function

## Details

This function firstly divides the individuals according to their
group/cluster membership, extracting group customized dataset and other
variables, and then compute the group log-likelihood function through
the function @ll_AdPaik_centre_1D. The produced group log-likelihood
value is summed together the other values into a unique result, that
corresponds to the overall (and final) log-likelihood value.
