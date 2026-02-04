# Plot the Posterior Frailty Estimates

This function plots the posterior frailty estimates for each group in
each time interval (represented by its mid point). Each group's
estimates are represented by a sequence of points connected by straight
lines. The function can plot either the entire posterior frailty
estimate or its time-independent and time-dependent components based on
user-specified flags.

## Usage

``` r
plot_post_frailty_est(
  result,
  flag_eps = FALSE,
  flag_alpha = FALSE,
  xlim = NULL,
  ylim = NULL,
  xlab = "Time",
  ylab = "Values",
  main = "Posterior frailty estimates",
  cex = 0.7,
  pch_type = seq(1, length(result$ClusterCodes)),
  color_bg = rep("black", length(result$ClusterCodes)),
  cex_legend = 0.7,
  pos_legend = "topright"
)
```

## Arguments

- result:

  S3 object of class 'AdPaik', returned by the method call
  'AdPaikModel(...)'.

- flag_eps:

  Logical flag indicating whether to plot only the time-dependent
  posterior frailty estimates. Default is FALSE.

- flag_alpha:

  Logical flag indicating whether to plot only the time-independent
  posterior frailty estimates. Default is FALSE.

- xlim:

  A numeric vector specifying the range for the x-axis (intervals). If
  NULL, default is set to the interval min-max of the time-domain, plus
  space for the legend. If flag_alpha = TRUE, the plot is produced
  around 1 (defaults to 0.8-1.4).

- ylim:

  A numeric vector specifying the range for the y-axis (intervals). If
  NULL, default is min-max value of the posterior frailty estimate.

- xlab, ylab:

  String giving the x and y axis name. Default values are 'Time' and
  'Values'.

- main:

  Title of the plot. Default title is 'Posterior frailty estimates'.

- cex:

  Dimension of the points used for plotting the estimates.

- pch_type:

  Numerical vector of length equal to the number of clusters in the
  data, giving the symbol to be used for plotting the estimates. Default
  symbol (circle, 21) is the same for all clusters.

- color_bg:

  Numerical vector of length equal to the number of clusters in the
  data, giving the color to be used for plotting the symbols for the
  estimates. Default ('black') is the same for all faculties. On the
  other hand, the same color is used throughout the intervals for the
  same faculty.

- cex_legend:

  Dimension of the symbol in the legend. Default is 0.7.

- pos_legend:

  Either a numeric vector providing the x and y coordinates for the
  legend or a string specifying the legend's position (e.g.,
  'bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top',
  'topright', 'right', 'center').

## Value

The plot of the posterior frailty estimates.

## Details

Recalling the frailty structure as \\Z\_{jk} = \alpha\_{j} +
\epsilon\_{jk}, \forall j,k\\ and the posterior frailty estimate as
\\\hat{Z}\_{jk} = \hat{\alpha}\_{j}/\hat{\alpha}\_{max} +
\hat{\epsilon}\_{jk}/\hat{\epsilon}\_{max}\\, this function allows
plotting either the entire posterior frailty estimate \\\hat{Z}\_{jk}\\
or its time-independent
\\\frac{\hat{\alpha}\_{j}}{\hat{\alpha}\_{\text{max}}}\\ or
time-dependent
\\\frac{\hat{\epsilon}\_{jk}}{\hat{\epsilon}\_{\text{max}}}\\
components. The user can control which components to display using the
flag_eps and flag_alpha parameters. Only one of these flags can be set
to TRUE at a time.

## Examples

``` r
# Import data
data(data_dropout)

# Define the variables needed for the model execution
eps_paik <- 1e-10
categories_range_min <- c(-8, -2, eps_paik, eps_paik, eps_paik)
categories_range_max <- c(-eps_paik, 0.4, 1 - eps_paik, 1, 10)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
formula <- time_to_event ~ Gender + CFUP + cluster(group)

# Call the main model function

# \donttest{
result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#> Error in while (r <= n_run & actual_tol_ll > tol_ll) {    if (verbose)         message(paste("Run ", r))    RemainingIndexes <- RunIndexes[r, ]    UsedIndexes <- c()    while (length(RemainingIndexes) != 0) {        index_to_vary <- RemainingIndexes[1]        PosIndex <- which(RemainingIndexes == index_to_vary)        RemainingIndexes <- RemainingIndexes[-PosIndex]        UsedIndexes <- c(UsedIndexes, index_to_vary)        result_optimize <- suppressWarnings(optimize(ll_AdPaik_1D,             c(params_range_min[index_to_vary], params_range_max[index_to_vary]),             maximum = TRUE, tol = tol_optimize, index_to_vary,             params, dataset, centre, time_axis, dropout_matrix,             e_matrix))        params[index_to_vary] <- result_optimize$maximum    }    global_optimal_params[r, ] <- params    global_optimal_loglikelihood_run <- ll_AdPaik_eval(params,         dataset, centre, time_axis, dropout_matrix, e_matrix)    global_optimal_loglikelihood[r] <- global_optimal_loglikelihood_run    if (is.nan(global_optimal_loglikelihood_run))         stop("NaN value for the optimal log-likelihood value.")    if (print_previous_ll_values[1]) {        n_previous <- print_previous_ll_values[2]        if (r < n_previous)             if (verbose)                 message(paste(" Global log-likelihood: ", global_optimal_loglikelihood[1:r]))            else if (verbose)                 message(paste(" Global log-likelihood: ", global_optimal_loglikelihood[(r -                   n_previous + 1):r]))    }    actual_tol_ll <- abs(ll_optimal - global_optimal_loglikelihood_run)    if (ll_optimal < global_optimal_loglikelihood_run) {        ll_optimal <- global_optimal_loglikelihood_run        optimal_run <- r    }    r <- r + 1}: missing value where TRUE/FALSE needed

# Define variables for plotting the estimates
pch_type <- c(21, seq(21,25,1), seq(21,25,1), seq(21,25,1))
color_bg <- c("darkblue", rep("red", 5), rep("purple", 5), rep("green",5))

plot_post_frailty_est(result, pch_type = pch_type, color_bg = color_bg)
#> Error: object 'result' not found
 # }                     
```
