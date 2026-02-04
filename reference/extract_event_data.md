# Extracting Variables for Posterior Frailty Estimates Computation

Function for extracting from the dataset quantities necessary to the
evaluation of the posterior frailty estimates.

## Usage

``` r
extract_event_data(dataset, time_to_event, centre, time_axis, phi, betar)
```

## Arguments

- dataset:

  Dataset containing the covariates/regressors. Their numerosity is
  indicated with R.

- time_to_event:

  Time-instant in the follow-up in which an individual fails or faces
  the event. If an individual does not face the event, the time-instant
  assumes a default value.

- centre:

  Categorical vector indicating the group/cluster membership. The number
  of distinct group is indicated with N.

- time_axis:

  Numerical vector of the temporal domain. Its length is (L+1), where L
  indicates the number of intervals of the time-domain.

- phi:

  Numerical vector of length L, of estimated baseline log-hazard.

- betar:

  Numerical vector of length R, of estimated regressors.

## Value

S3 object of class 'EventData', composed of six elements. See details.

## Details

The S3 class obejct 'EventData' contains the variables necessary for the
estimate of the posterior frailty and that can be extracted or computed
starting from the dataset.

- N_ik: matrix of dimension (N, L), containing the number of event in
  each interval k and group i.

- N_i: numerical vector of length L, with the number of event in each
  group i. It can be computed as: \\N_i = \sum\_{k=1}^L N\_{ik}\\.

- e_ijk: matrix of dimension (n_individuals, L) with the evaluation of
  the temporal integral, for each individual j, group i and interval k.

- Y_risk: binary matrix of dimension (n_individuals, L) reporting for
  each individual, in each interval, his/her risk of facing the event.
  For an individual, the risk is equal to 1 in an interval k if, in that
  interval, he/she has not faced the event yet; otherwise, it is equal
  to 0.

- cum_hazard_group: matrix of dimension (N, L), where each element in
  position (i,k) indicates the computed cumulative hazard for all
  individuals belonging to group i and at interval k.

- sum_cum_hazard_group: numerical vector of length N, giving the sum of
  the computed cumulative hazard for all intervals k and for all
  individuals belonging to group i. It can be computed from the previous
  element, summing with respect to the interval k.
