# Frailty Standard Deviation and Variance for the 'Adapted Paik et Al.'s Model'

The function computes both the standard deviation and the variance of
the time-dependent frailty of the 'Adapted Paik et al.'s Model'.

Recalling the frailty structure \\Z\_{jk} = \alpha_j + \epsilon\_{jk}\\
as being composed by a constant group-dependent term (\\\alpha_j\\) and
a time and group dependent term (\\\epsilon\_{jk}\\), the frailty
variance (and standard deviation) can be computed in two different way:

- Considering only the time-dependent spread of the
  clusters/groups/centre: \\var(Z\_{jk}) = \mu_2 \* \gamma_k\\. In this
  case, the flag_full should be FALSE and flag_variance should be TRUE.

- Considering both the time-dependent and constant spread of the
  clusters: \\var(Z\_{jk}) = \mu_1 \* \nu + \mu_2 \* \gamma_k\\. The new
  added term only moves upward the other case and the flag_full should
  be TRUE and flag_variance should be TRUE.

## Usage

``` r
frailty_sd(object, flag_full = TRUE, flag_variance = FALSE)
```

## Arguments

- object:

  S3 object of class 'AdPaik' returned by the main model output, that
  contains all the information for the computation of the frailty
  standard deviation.

- flag_full:

  A boolean flag indicating whether to get the full standard deviation
  (`TRUE`) or only the time-dependent component (`FALSE`). Default to
  `TRUE`.

- flag_variance:

  A boolean flag indicating whether to get the frailty variance (`TRUE`)
  or the frailty standard deviation (`FALSE`). Default to `FALSE`.

## Value

Numerical vector of length equal to the number of intervals of the
time-domain, with the value of the frailty standard deviation or
variance (either full or only the time-dependent component).

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

frailty_sd(result)
#>  [1] 0.07866498 0.24103472 0.27738443 0.07701667 0.12557314 0.28758349
#>  [7] 0.42205690 0.42085975 0.07678833 0.10058504
# }
```
