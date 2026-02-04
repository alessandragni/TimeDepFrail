# Plot the Baseline Hazard Step-Function

This function plots the baseline hazard step-function based on the
estimated parameters from the Adapted Paik et al.'s model.

## Usage

``` r
plot_bas_hazard(
  result,
  xlim = c(min(result$TimeDomain), max(result$TimeDomain)),
  ylim = c(0, max(result$BaselineHazard)),
  xlab = "Time",
  ylab = "Values",
  main = "Baseline hazard step-function",
  color = "black",
  pch = 21,
  bg = "black",
  cex_points = 0.7
)
```

## Arguments

- result:

  S3 object of class 'AdPaik', returned by the method call
  'AdPaikModel(...)'.

- xlim:

  A numeric vector specifying the x-axis limits. Default is set to the
  interval min-max of the time-domain.

- ylim:

  A numeric vector specifying the y-axis limits. Default is 0 to the
  maximum value of the baseline hazard.

- xlab, ylab:

  String giving the x and y axis name. Default values are 'x' and 'y'.

- main:

  Title of the plot. Default title is 'Baseline hazard step-function'.

- color:

  Color used for plotting the horizontal segments of the step-function.
  Default one is 'black'.

- pch:

  Symbol for marking the boundaries of each segment. Default is a dot
  (value 21).

- bg:

  Color for the boundary symbols. Default matches the plot color
  ('black').

- cex_points:

  Size of the boundary symbols. Default is 0.7.

## Value

Plot of the baseline hazard step-function and value of the function in
each interval.

## Details

The function plots a horizontal segment for each interval of the time
domain, representing the baseline hazard. The boundaries of each segment
are marked with colored dots, and subsequent segments are intentionally
left unconnected to reflect the discrete nature of the intervals.

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

# \donttest{
# Call the main model function
result <- AdPaikModel(formula, data_dropout, time_axis, categories_range_min, categories_range_max)
#> Error in while (r <= n_run & actual_tol_ll > tol_ll) {    if (verbose)         message(paste("Run ", r))    RemainingIndexes <- RunIndexes[r, ]    UsedIndexes <- c()    while (length(RemainingIndexes) != 0) {        index_to_vary <- RemainingIndexes[1]        PosIndex <- which(RemainingIndexes == index_to_vary)        RemainingIndexes <- RemainingIndexes[-PosIndex]        UsedIndexes <- c(UsedIndexes, index_to_vary)        result_optimize <- suppressWarnings(optimize(ll_AdPaik_1D,             c(params_range_min[index_to_vary], params_range_max[index_to_vary]),             maximum = TRUE, tol = tol_optimize, index_to_vary,             params, dataset, centre, time_axis, dropout_matrix,             e_matrix))        params[index_to_vary] <- result_optimize$maximum    }    global_optimal_params[r, ] <- params    global_optimal_loglikelihood_run <- ll_AdPaik_eval(params,         dataset, centre, time_axis, dropout_matrix, e_matrix)    global_optimal_loglikelihood[r] <- global_optimal_loglikelihood_run    if (is.nan(global_optimal_loglikelihood_run))         stop("NaN value for the optimal log-likelihood value.")    if (print_previous_ll_values[1]) {        n_previous <- print_previous_ll_values[2]        if (r < n_previous)             if (verbose)                 message(paste(" Global log-likelihood: ", global_optimal_loglikelihood[1:r]))            else if (verbose)                 message(paste(" Global log-likelihood: ", global_optimal_loglikelihood[(r -                   n_previous + 1):r]))    }    actual_tol_ll <- abs(ll_optimal - global_optimal_loglikelihood_run)    if (ll_optimal < global_optimal_loglikelihood_run) {        ll_optimal <- global_optimal_loglikelihood_run        optimal_run <- r    }    r <- r + 1}: missing value where TRUE/FALSE needed

plot_bas_hazard(result)
#> Error: object 'result' not found
# }
```
