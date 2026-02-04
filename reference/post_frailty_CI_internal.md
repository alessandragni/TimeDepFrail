# Confidence Interval for Posterior Frailty Estimates

Function for computing the confidence interval for each posterior
frailty estimates \\\hat{Z}\_{jk}\\.

## Usage

``` r
post_frailty_CI_internal(
  post_frailty_est,
  post_frailty_est_var,
  n_centres,
  n_intervals,
  level
)
```

## Arguments

- post_frailty_est:

  Posterior frailty estimates list.

- post_frailty_est_var:

  Posterior frailty variance list.

- n_centres:

  Number of clusters/centres.

- n_intervals:

  Number of intervals of the time-domain. it is equal to the length of
  the tima_axis minus one.

- level:

  A numeric value representing the confidence level.

## Value

S3 object of class 'PFCI.AdPaik' composed of two matrices of dimension
(number groups, number of intervals):

- PostFrailtyCI_left: left confidence interval for each posterior
  frailty estimates

- PostFrailtyCI_right: right confidence interval for each each posterior
  frailty estimates
