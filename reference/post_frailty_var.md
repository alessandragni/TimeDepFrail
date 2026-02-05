# Posterior Frailty Variances

Function for computing the posterior frailty variances of the
time-dependent shared frailty Cox model.

Recalling the structure of the frailty \\Z\_{jk} = \alpha_j +
\epsilon\_{jk}, \forall j,k\\ with \\k=1,\dots,L\\ and \\j=1,\dots,N\\
as being composed by the sum of two independent gamma distributions:

- \\\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j\\

- \\\epsilon\_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k\\

The posterior frailty variance is \\var(\hat{Z}\_{jk}) =
var(\hat{\alpha}\_{j}/\hat{\alpha}\_{max}) +
var(\hat{\epsilon}\_{jk}/\hat{\epsilon}\_{max}\\). This function allows
to get either the entire posterior frailty variance
\\var(\hat{Z}\_{jk})\\ or its time-independent
\\var(\frac{\hat{\alpha}\_{j}}{\hat{\alpha}\_{\text{max}}})\\ or
time-dependent
\\var(\frac{\hat{\epsilon}\_{jk}}{\hat{\epsilon}\_{\text{max}}})\\
components. The user can control which components to display using the
flag_eps and flag_alpha parameters. Only one of these flags can be set
to TRUE at a time.

## Usage

``` r
post_frailty_var(object, flag_eps = FALSE, flag_alpha = FALSE)
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

Vector or matrix of posterior frailty variances, depending on the
flag_eps and flag_alpha values. Specifically:

- It is a vector of length equal to the N containing posterior frailty
  variances for \\\alpha_j, \forall j\\. In this case the flag_eps must
  be FALSE and the flag_alpha must be TRUE.

- Matrix of dimension (N, L) containing posterior frailty variances for
  \\\epsilon\_{jk}, \forall j,k\\. In this case the flag_eps must be
  TRUE and the flag_alpha must be FALSE.

- Matrix of dimension (N, L) containing posterior frailty variances for
  \\Z\_{jk} \forall j,k\\. In this case the flag_eps must be FALSE and
  the flag_alpha must be FALSE.

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

post_frailty_var(result)
#>              [,1]       [,2]       [,3]        [,4]        [,5]       [,6]
#>  [1,] 0.009208063 0.01948519 0.01625035 0.009189044 0.011345676 0.01675327
#>  [2,] 0.008539039 0.01817861 0.01315513 0.008519156 0.010796773 0.01997050
#>  [3,] 0.035183632 0.04655819 0.04659670 0.035128557 0.037541483 0.05106701
#>  [4,] 0.016448517 0.03043441 0.02706092 0.016410353 0.018913295 0.02753686
#>  [5,] 0.012632935 0.02382874 0.01908481 0.012583300 0.014817080 0.02400076
#>  [6,] 0.016924983 0.02747806 0.02728709 0.016921081 0.019342300 0.02653674
#>  [7,] 0.046995425 0.05971937 0.05902819 0.046932212 0.049192633 0.06076348
#>  [8,] 0.022012324 0.03307887 0.03135589 0.021937808 0.024141574 0.03523329
#>  [9,] 0.032894965 0.04480348 0.04406633 0.032844002 0.035310595 0.04788441
#> [10,] 0.008963879 0.01952475 0.01473208 0.008930203 0.010906160 0.01790961
#> [11,] 0.033919786 0.04516761 0.04968049 0.033862364 0.036117181 0.04895887
#> [12,] 0.008987804 0.01912232 0.01576456 0.009031575 0.011211195 0.02222004
#> [13,] 0.012685761 0.02151768 0.02026928 0.012672314 0.014739978 0.02456902
#> [14,] 0.007730287 0.01592580 0.01184410 0.007712729 0.009728632 0.01914208
#> [15,] 0.025736224 0.03710276 0.03615343 0.025683808 0.027970719 0.04124599
#> [16,] 0.016810834 0.02804933 0.02583542 0.016752506 0.019091456 0.03049263
#>             [,7]       [,8]        [,9]       [,10]
#>  [1,] 0.03837373 0.02367612 0.009189142 0.010243468
#>  [2,] 0.04191872 0.02535611 0.008509982 0.009578139
#>  [3,] 0.07379936 0.05902863 0.035117215 0.036139735
#>  [4,] 0.04807015 0.03201288 0.016422521 0.017510260
#>  [5,] 0.04595105 0.03762439 0.012581302 0.013624871
#>  [6,] 0.05096962 0.03837089 0.016909497 0.017971329
#>  [7,] 0.08492578 0.08396340 0.046904799 0.048011215
#>  [8,] 0.06461598 0.06305206 0.021941110 0.022993019
#>  [9,] 0.07087392 0.06905323 0.032843144 0.033924159
#> [10,] 0.05594357 0.02253692 0.008946395 0.010009140
#> [11,] 0.07109628 0.06700955 0.033871019 0.034916399
#> [12,] 0.04395489 0.02000222 0.009037689 0.010126022
#> [13,] 0.05655574 0.05180134 0.012670481 0.013752401
#> [14,] 0.04146032 0.01866148 0.007745280 0.008937592
#> [15,] 0.06928247 0.06061712 0.025676401 0.026703864
#> [16,] 0.05141540 0.03347129 0.016757913 0.017853681
# }
```
