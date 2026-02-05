# Plot for the Frailty Standard Deviation or Variance

This function generates a plot of either the frailty standard deviation
or the frailty variance for the intervals in the time-domain.

## Usage

``` r
plot_frailty_sd(
  result,
  flag_full = TRUE,
  flag_variance = FALSE,
  xlim = c(min(result$TimeDomain), max(result$TimeDomain)),
  ylim = NULL,
  xlab = "Time",
  ylab = "Values",
  main = NULL,
  pch = 21,
  color_bg = "blue",
  cex_points = 0.7
)
```

## Arguments

- result:

  An S3 object of class 'AdPaik', returned by the main model call
  'AdPaikModel(...)'.

- flag_full:

  A boolean flag indicating whether to plot the full standard deviation
  (`TRUE`) or only the time-dependent one (`FALSE`). Default is `TRUE`.

- flag_variance:

  A boolean flag indicating whether to plot the frailty variance
  (`TRUE`) or the frailty standard deviation (`FALSE`). Default is
  `FALSE`.

- xlim:

  A numeric vector specifying the range for the x-axis (intervals). If
  NULL, default is set to the interval min-max of the time-domain.

- ylim:

  A numeric vector specifying the range for the y-axis (intervals). If
  NULL, default is 0 to the maximum value of the frailty
  variance/standard deviation.

- xlab:

  A string for the x-axis label. Default is `'Intervals'`.

- ylab:

  A string for the y-axis label. Default is `'Values'`.

- main:

  A string for the plot title. Default title is
  `'Frailty Standard Deviation'` or `'Frailty Variance'` according to
  the produced plot (flag_variance).

- pch:

  A numeric or character symbol used for plotting the frailty standard
  deviation values. Default is a dot (`21`).

- color_bg:

  A string specifying the color used for the symbols. Default is
  `'blue'`.

- cex_points:

  A numeric value indicating the size of the plotting symbols. Default
  is `0.7`.

## Value

A plot displaying either the frailty standard deviation or variance
across the specified intervals.

## Details

The plot represents the values of the frailty standard deviation or
variance for each time interval (represented by its mid point). It
connects these points to illustrate the trend of the chosen metric.

This function supports plotting the full or only time dependent frailty
standard deviation or variance retrieved from the main model (contained
in the S3 object of class 'AdPaik').

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

plot_frailty_sd(result)

# }
```
