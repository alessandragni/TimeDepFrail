# Check Coherence Between Flag for Optimal Parameters and Optimal Parameters

The function controls that one of the following condition is satisfied:

- if the flag for the optimal parameters is activated, then the optimal
  parameters should be provided in input

- if the flag is not activated, then the optimal parameters should not
  be provided and the parameter vector should be NU

## Usage

``` r
check.flag_optimal_params(optimal_params, flag_optimal_params)
```

## Arguments

- optimal_params:

  Either a numerical vector of length equal to the number of model
  parameters or a NULL value.

- flag_optimal_params:

  Logical value. Did the user want to provide optimal parameters vector?
  If so, the variable should be TRUE; otherwise (no optimal parameters),
  FALSE.

## Value

An error if any condition is not satisfied.
