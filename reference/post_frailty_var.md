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
#>  [1,] 0.009208125 0.01948534 0.01625038 0.009189107 0.011345461 0.01675335
#>  [2,] 0.008539112 0.01817878 0.01315520 0.008519229 0.010796548 0.01997065
#>  [3,] 0.035183624 0.04655832 0.04659671 0.035128549 0.037541154 0.05106715
#>  [4,] 0.016448597 0.03043463 0.02706096 0.016410434 0.018913045 0.02753699
#>  [5,] 0.012633024 0.02382895 0.01908489 0.012583388 0.014816879 0.02400093
#>  [6,] 0.016924958 0.02747814 0.02728701 0.016921057 0.019341952 0.02653675
#>  [7,] 0.046995443 0.05971954 0.05902822 0.046932230 0.049192353 0.06076362
#>  [8,] 0.022012527 0.03307920 0.03135611 0.021938010 0.024141493 0.03523361
#>  [9,] 0.032895158 0.04480381 0.04406655 0.032844195 0.035310461 0.04788475
#> [10,] 0.008963939 0.01952489 0.01473210 0.008930264 0.010905970 0.01790969
#> [11,] 0.033920030 0.04516799 0.04968077 0.033862608 0.036117134 0.04895926
#> [12,] 0.008987855 0.01912244 0.01576455 0.009031628 0.011210966 0.02222010
#> [13,] 0.012685842 0.02151783 0.02026932 0.012672396 0.014739793 0.02456915
#> [14,] 0.007730331 0.01592589 0.01184411 0.007712774 0.009728421 0.01914213
#> [15,] 0.025736333 0.03710299 0.03615355 0.025683916 0.027970529 0.04124623
#> [16,] 0.016810925 0.02804954 0.02583550 0.016752597 0.019091241 0.03049282
#>             [,7]       [,8]        [,9]       [,10]
#>  [1,] 0.03837415 0.02367618 0.009189205 0.010243553
#>  [2,] 0.04191927 0.02535625 0.008510056 0.009578235
#>  [3,] 0.07379996 0.05902874 0.035117206 0.036139749
#>  [4,] 0.04807065 0.03201297 0.016422602 0.017510364
#>  [5,] 0.04595162 0.03762457 0.012581390 0.013624982
#>  [6,] 0.05097005 0.03837086 0.016909473 0.017971328
#>  [7,] 0.08492640 0.08396360 0.046904816 0.048011256
#>  [8,] 0.06461688 0.06305252 0.021941312 0.022993244
#>  [9,] 0.07087472 0.06905367 0.032843337 0.033924375
#> [10,] 0.05594427 0.02253696 0.008946456 0.010009224
#> [11,] 0.07109711 0.06701000 0.033871264 0.034916666
#> [12,] 0.04395535 0.02000224 0.009037742 0.010126099
#> [13,] 0.05655643 0.05180140 0.012670563 0.013752506
#> [14,] 0.04146075 0.01866149 0.007745326 0.008937664
#> [15,] 0.06928329 0.06061740 0.025676509 0.026703995
#> [16,] 0.05141600 0.03347143 0.016758004 0.017853796
# }
```
