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
#> Error in ll_AdPaik_centre_1D(x, index, params, dataset_centre, dropout_matrix_centre,     e_matrix_centre): dims [product 330] do not match the length of object [3300]

post_frailty_var(result)
#> Error: object 'result' not found
# }
```
