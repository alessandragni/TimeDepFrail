# Posterior Frailty Confidence Intervals

Function for computing the posterior frailty confidence intervals of the
time-dependent shared frailty Cox model.

Recalling the structure of the frailty \\Z\_{jk} = \alpha_j +
\epsilon\_{jk}, \forall j,k\\ with \\k=1,\dots,L\\ and \\j=1,\dots,N\\
as being composed by the sum of two independent gamma distributions:

- \\\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j\\

- \\\epsilon\_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k\\

The posterior frailty estimate is \\\hat{Z}\_{jk} =
\hat{\alpha}\_{j}/\hat{\alpha}\_{max} +
\hat{\epsilon}\_{jk}/\hat{\epsilon}\_{max}\\. This function allows to
get the confidence intervals of either the entire posterior frailty
estimates \\\hat{Z}\_{jk}\\ or its time-independent
\\\frac{\hat{\alpha}\_{j}}{\hat{\alpha}\_{\text{max}}}\\ or
time-dependent
\\\frac{\hat{\epsilon}\_{jk}}{\hat{\epsilon}\_{\text{max}}}\\
components. The user can control which components to display using the
flag_eps and flag_alpha parameters. Only one of these flags can be set
to TRUE at a time.

## Usage

``` r
post_frailty_confint(
  object,
  level = 0.95,
  flag_eps = FALSE,
  flag_alpha = FALSE
)
```

## Arguments

- object:

  S3 object of class 'AdPaik' returned by the main model output, that
  contains all the information for the computation of the frailty
  standard deviation.

- level:

  A numeric value representing the confidence level for the posterior
  frailty confidence intervals. Default is 0.95 for 95% confidence.

- flag_eps:

  Logical flag indicating whether to extract only the time-dependent
  posterior frailty estimates. Default is FALSE.

- flag_alpha:

  Logical flag indicating whether to extract only the time-independent
  posterior frailty estimates. Default is FALSE.

## Value

A list for posterior frailty confidence intervals, depending on the
flag_eps and flag_alpha values. Specifically:

- A list of length equal to the N containing posterior frailty
  confidence intervals for \\\alpha_j, \forall j\\. In this case the
  flag_eps must be FALSE and the flag_alpha must be TRUE.

- A list of length equal to the NxL containing posterior frailty
  confidence intervals for \\\epsilon\_{jk}, \forall j,k\\. In this case
  the flag_eps must be TRUE and the flag_alpha must be FALSE.

- A list of length equal to the NxL containing posterior frailty
  confidence intervals for \\Z\_{jk} \forall j,k\\. In this case the
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
#> Error in while (r <= n_run & actual_tol_ll > tol_ll) {    if (verbose)         message(paste("Run ", r))    RemainingIndexes <- RunIndexes[r, ]    UsedIndexes <- c()    while (length(RemainingIndexes) != 0) {        index_to_vary <- RemainingIndexes[1]        PosIndex <- which(RemainingIndexes == index_to_vary)        RemainingIndexes <- RemainingIndexes[-PosIndex]        UsedIndexes <- c(UsedIndexes, index_to_vary)        result_optimize <- suppressWarnings(optimize(ll_AdPaik_1D,             c(params_range_min[index_to_vary], params_range_max[index_to_vary]),             maximum = TRUE, tol = tol_optimize, index_to_vary,             params, dataset, centre, time_axis, dropout_matrix,             e_matrix))        params[index_to_vary] <- result_optimize$maximum    }    global_optimal_params[r, ] <- params    global_optimal_loglikelihood_run <- ll_AdPaik_eval(params,         dataset, centre, time_axis, dropout_matrix, e_matrix)    global_optimal_loglikelihood[r] <- global_optimal_loglikelihood_run    if (is.nan(global_optimal_loglikelihood_run))         stop("NaN value for the optimal log-likelihood value.")    if (print_previous_ll_values[1]) {        n_previous <- print_previous_ll_values[2]        if (r < n_previous)             if (verbose)                 message(paste(" Global log-likelihood: ", global_optimal_loglikelihood[1:r]))            else if (verbose)                 message(paste(" Global log-likelihood: ", global_optimal_loglikelihood[(r -                   n_previous + 1):r]))    }    actual_tol_ll <- abs(ll_optimal - global_optimal_loglikelihood_run)    if (ll_optimal < global_optimal_loglikelihood_run) {        ll_optimal <- global_optimal_loglikelihood_run        optimal_run <- r    }    r <- r + 1}: missing value where TRUE/FALSE needed

post_frailty_confint(result)
#> Error: object 'result' not found
# }
```
