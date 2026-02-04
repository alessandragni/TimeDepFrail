# Check Correctness of Plot Variables Pch and Color

The function controls that the input variables 'pch_type' and 'color_bg'
have the correct structure, they have the same dimension of the number
of clusters in the dataset and they have meaningful elements.

These variables are used for the plot of the posterio frailty estimates:
the estimates for each faculty are plotted through a symbol, having
color and shape indicated by the variables (for the k-th faculty,
consider the k-th element of both vectors).

## Usage

``` r
check.pchtype_colorbg(centre_codes, pch_type, color_bg)
```

## Arguments

- centre_codes:

  Numerical vector of length equal to the number of centres/clusters in
  the dataset and containing the distinct centres/clusters. They
  correspond to the levels of the numerical vector of individual group
  membership.

- pch_type:

  Numerical vector of length equal to the number of centres and
  containing the point shape for each faculty.

- color_bg:

  Numerical vector of length equal to the number of centres and
  containing the color of the point for each faculty.

## Value

An error if any condition is not satisfied.
