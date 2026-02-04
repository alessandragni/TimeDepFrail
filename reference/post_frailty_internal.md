# Posterior Frailty Estimates and Variances for the 'Adapted Paik et Al.'s Model'

Posterior Frailty Estimates and Variances for the 'Adapted Paik et Al.'s
Model'

## Usage

``` r
post_frailty_internal(
  optimal_params,
  dataset,
  time_to_event,
  centre,
  time_axis
)
```

## Arguments

- optimal_params:

  Optimal parameters estimated by maximizing the log-likelihood
  function, through the constraint multi-dimensional optmization method.

- dataset:

  Dataset containing all the covariates/regressors.

- time_to_event:

  Time-instant, in the follow-up, in which an individual faces the event
  or fails. If an individual does not face the event in the follow-up,
  then the time-instant must assume a default value.

- centre:

  Individual group/cluster membership.

- time_axis:

  Temporal domain.

## Value

S3 object of class 'PF.AdPaik' composed of two elements of different
class:

- PosteriorFrailtyEst: S3 object of class 'PFE.AdPaik'.

- PosteriorFrailtyVar: S3 object of class 'PFV.AdPaik'.
