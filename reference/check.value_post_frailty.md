# Check Non-Negativeness of the Posterior Frailty Estimates

The function controls that all posterior frailty estimates are
non-negative. Indeed, by construction the realizations of a gamma
distribution are non negative.

## Usage

``` r
check.value_post_frailty(post_frailty_est, n_centres, n_intervals)
```

## Arguments

- post_frailty_est:

  An S3 class object containing the posterior frailty estimates:
  \\\hat{\alpha}\_j, \hat{\epsilon}\_{jk}, \hat{Z}\_{jk}, \forall j,k\\

- n_centres:

  Number of groups/clusters.

- n_intervals:

  Number of intervals of the time domain

## Value

An error if any condition is not satisfied.
