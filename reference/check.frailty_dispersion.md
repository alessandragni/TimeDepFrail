# Check Correctness of Frailty Standard Deviation

The function controls that the frailty standard deviation vector has a
length equal to the number of inyervals of the time domain and that its
elements are non-negative.

## Usage

``` r
check.frailty_dispersion(frailty_dispersion, n_intervals)
```

## Arguments

- frailty_dispersion:

  Frailty dispersion

- n_intervals:

  Number of intervals of the time-domain

## Value

An error if any condition is not satisfied.
