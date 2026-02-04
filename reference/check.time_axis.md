# Check Correctness of Time Domain Subdivision

The function controls that the time domain is a vector and it has at
least 2 elements and that all of them are not negative. Moreover, it
checks that all its elements are non-negative and properly ordered, in
an ascending way.

## Usage

``` r
check.time_axis(time_axis)
```

## Arguments

- time_axis:

  Numerical vector of temporal domain subdivision.

## Value

An error is returned if any condition is not satisfied.
