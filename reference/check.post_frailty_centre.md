# Check Numerosity of Posterior Frailty Estimates

The function controls that a time-dependent posterior frailty estimate
is computed for each centre

## Usage

``` r
check.post_frailty_centre(post_frailty_est, centre_codes)
```

## Arguments

- post_frailty_est:

  An S3 class object containing the posterior frailty estimates
  \\\hat{\alpha}\_j, \hat{\epsilon}\_{jk}, \hat{Z}\_{jk}, \forall j,k\\

- centre_codes:

  Numerical vector of length equal to the number of distinct
  centres/clusters in the study

## Value

An error if any condition is not satisfied.
