# Posterior Frailty Estimates

Function for computing the posterior frailty estimates of the
time-dependent shared frailty Cox model.

Recalling the structure of the frailty \\Z\_{jk} = \alpha_j +
\epsilon\_{jk}, \forall j,k\\ with \\k=1,\dots,L\\ and \\j=1,\dots,N\\
as being composed by the sum of two independent gamma distributions:

- \\\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j\\

- \\\epsilon\_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k\\

The posterior frailty estimate is \\\hat{Z}\_{jk} =
\hat{\alpha}\_{j}/\hat{\alpha}\_{max} +
\hat{\epsilon}\_{jk}/\hat{\epsilon}\_{max}\\. This function allows to
get either the entire posterior frailty estimate \\\hat{Z}\_{jk}\\ or
its time-independent
\\\frac{\hat{\alpha}\_{j}}{\hat{\alpha}\_{\text{max}}}\\ or
time-dependent
\\\frac{\hat{\epsilon}\_{jk}}{\hat{\epsilon}\_{\text{max}}}\\
components. The user can control which components to display using the
flag_eps and flag_alpha parameters. Only one of these flags can be set
to TRUE at a time.

## Usage

``` r
post_frailty_est(object, flag_eps = FALSE, flag_alpha = FALSE)
```

## Arguments

- object:

  S3 object of class 'AdPaik' returned by the main model output, that
  contains all the information for the computation of the frailty
  standard deviation.

- flag_eps:

  Logical flag indicating whether to extract only the time-dependent
  posterior frailty estimates. Default is FALSE.

- flag_alpha:

  Logical flag indicating whether to extract only the time-independent
  posterior frailty estimates. Default is FALSE.

## Value

Vector or matrix of posterior frailty estimates, depending on the
flag_eps and flag_alpha values. Specifically:

- It is a vector of length equal to the N containing posterior frailty
  estimates for \\\alpha_j, \forall j\\. In this case the flag_eps must
  be FALSE and the flag_alpha must be TRUE.

- Matrix of dimension (N, L) containing posterior frailty estimates for
  \\\epsilon\_{jk}, \forall j,k\\. In this case the flag_eps must be
  TRUE and the flag_alpha must be FALSE.

- Matrix of dimension (N, L) containing posterior frailty estimates for
  \\Z\_{jk} = \alpha_j + \epsilon\_{jk}, \forall j,k\\. In this case the
  flag_eps must be FALSE and the flag_alpha must be FALSE.

## Examples

``` r
# Consider the 'Academic Dropout dataset'
data(data_dropout)

# Define the variables needed for the model execution
formula <- time_to_event ~ Gender + CFUP + cluster(group)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)

# \donttest{
# Call the main model
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)

post_frailty_est(result)
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 1.1048039 1.1087732 1.1474586 1.1059445 1.1023892 0.9904662 1.0244111
#>  [2,] 0.8378979 0.7839215 0.6487738 0.8448966 0.8444854 0.8022525 0.7888660
#>  [3,] 1.1625652 1.1326906 1.1164246 1.1623098 1.1654944 1.1567702 1.1361859
#>  [4,] 1.4271670 1.5685548 1.6006521 1.4236172 1.4622315 1.4087016 1.3630419
#>  [5,] 1.0151196 1.0107879 0.8951994 1.0126693 1.0077689 0.9719493 0.9590378
#>  [6,] 1.1586229 1.1374121 1.2133199 1.1713000 1.1852619 1.0503538 1.1159607
#>  [7,] 1.5066852 1.5235085 1.4942492 1.5035623 1.4868891 1.4476758 1.4707267
#>  [8,] 1.1898357 1.1653092 1.1408120 1.1802733 1.1676181 1.1581835 1.2267306
#>  [9,] 1.2029260 1.2029242 1.1728395 1.2033384 1.2188186 1.1907841 1.1795056
#> [10,] 1.2415504 1.2794234 1.2640982 1.2343000 1.2156332 1.2301228 1.4222327
#> [11,] 1.3800599 1.3639587 1.5389786 1.3763672 1.3666791 1.3816889 1.3479933
#> [12,] 1.3800353 1.4470644 1.5676158 1.3933804 1.4089559 1.6238949 1.4235839
#> [13,] 1.4914496 1.4500259 1.6247025 1.4920520 1.4806392 1.5839651 1.6188102
#> [14,] 1.2537506 1.2086601 1.2219552 1.2492583 1.2451371 1.4184806 1.2740265
#> [15,] 1.2117142 1.1975282 1.1892655 1.2103528 1.2043611 1.2365914 1.2591603
#> [16,] 1.1414413 1.1376630 1.1241080 1.1360517 1.1440767 1.1523408 1.0931520
#>            [,8]      [,9]    [,10]
#>  [1,] 1.0215688 1.1045113 1.111616
#>  [2,] 0.7205009 0.8418721 0.850116
#>  [3,] 1.0399464 1.1603806 1.156733
#>  [4,] 1.3158036 1.4271675 1.439955
#>  [5,] 1.0479792 1.0123338 1.015600
#>  [6,] 1.1145690 1.1678633 1.173345
#>  [7,] 1.5770085 1.4959602 1.510278
#>  [8,] 1.3905406 1.1823069 1.186246
#>  [9,] 1.2793360 1.2046627 1.213948
#> [10,] 1.1758823 1.2370915 1.245612
#> [11,] 1.4289020 1.3807795 1.381963
#> [12,] 1.2935066 1.3915388 1.407312
#> [13,] 1.9939448 1.4888332 1.501144
#> [14,] 1.1674298 1.2563265 1.290922
#> [15,] 1.3057394 1.2090939 1.208076
#> [16,] 0.9903548 1.1384444 1.151889
# }
```
