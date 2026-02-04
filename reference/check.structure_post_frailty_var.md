# Check Structure of Posterior Frailty Variances

The function controls that the structure of the 'Posterior Frailty
Variances' coincides with the theoretical one.

## Usage

``` r
check.structure_post_frailty_var(post_frailty_var, n_intervals, n_centres)
```

## Arguments

- post_frailty_var:

  Posterior frailty variances S3 object of class 'PFV.AdPaik', composed
  of three elements:

  - 'alphaVar': posterior frailty variance for \\\alpha_j, \forall j\\.
    It is a vector of length equal to the number of centres.

  - 'epsVar': posterior frailty variance for \\\epsilon\_{jk}, \forall
    j,k\\. It is a matrix of dimension (n_centres, n_intervals).

  - 'ZVar': posterior frailty variance for \\Z\_{jk} = \alpha_j +
    \epsilon\_{jk}, \forall j,k\\. It is a matrix of dimension
    (n_centres, n_intervals)

- n_intervals:

  Number of intervals of the time-domain

- n_centres:

  Number of centres/clusters.

## Value

An error if any condition is not satisfied.
