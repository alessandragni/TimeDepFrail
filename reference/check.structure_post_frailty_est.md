# Check Structure of Posterior Frailty Estimates

The function controls that the structure of the 'Posterior Frailty
Estimates' coincides with the theoretical one.

## Usage

``` r
check.structure_post_frailty_est(post_frailty_est, n_intervals, n_centres)
```

## Arguments

- post_frailty_est:

  Posterior frailty estimates S3 object of class 'PFE.AdPaik', composed
  of three elements:

  - 'alpha': posterior frailty estimates for \\\alpha_j, \forall j\\. It
    is a vector of length equal to the number of centres.

  - 'eps': posterior frailty estimates for \\\epsilon\_{jk}, \forall
    j,k\\. It is a matrix of dimension (n_centres, n_intervals).

  - 'Z': posterior frailty estimates for \\Z\_{jk} = \alpha_j +
    \epsilon\_{jk}, \forall j,k\\. It is a matrix of dimension
    (n_centres, n_intervals)

- n_intervals:

  Number of intervals of the time-domain

- n_centres:

  Number of centres/clusters.

## Value

An error if any condition is not satisfied.
