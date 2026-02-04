# Check Correctness for the Cluster Variable

The function controls that the provided cluster variable is a vector,
with at least two levels. Indeed, it is not possible to apply the
Time-Dependent Shared Frailty Cox Model with no clusters.

## Usage

``` r
check.centre(centre)
```

## Arguments

- centre:

  Numerical vector of length equal to the number of individuals in the
  study, containing the individual grouo/cluster membership.

## Value

An error if any condition is not satisfied.
