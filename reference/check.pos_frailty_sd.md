# Check Positivity of the Frailty Standard Deviation

The method controls that the frailty standard deviation vector has
non-negative elements

## Usage

``` r
check.pos_frailty_sd(sd, n_intervals)
```

## Arguments

- sd:

  Numerical vector of length equal to the number of intervals,
  containing the frailty standard deviation

- n_intervals:

  Number of intervals of the time-domain

## Value

An error if any condition is not satisfied.
