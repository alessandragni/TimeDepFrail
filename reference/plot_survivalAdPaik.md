# Plot of Conditional Survival Function

Plots the conditional survival function based on the 'Adapted Paik et
al.' model's estimated coefficients and frailty effects, for each unit
in each time interval (represented by its mid point).

## Usage

``` r
plot_survivalAdPaik(
  result,
  lwd = 1,
  xlim = c(min(result$TimeDomain), max(result$TimeDomain)),
  ylim = c(0, 1),
  xlab = "Time",
  ylab = "Values",
  main = "Conditional Survival",
  cex = 0.2,
  cexlegend = 0.8
)
```

## Arguments

- result:

  S3 object of class 'AdPaik' containing model results.

- lwd:

  The line width of the plot. Default is 1.

- xlim:

  A numeric vector specifying the range for the x-axis (intervals).
  Default is min-max value of the time domain.

- ylim:

  A numeric vector specifying the range for the y-axis (intervals).
  Default is the range 0-1.

- xlab, ylab:

  String giving the x and y axis name. Default values are 'Time' and
  'Values'.

- main:

  Title of the plot. Default title is 'Survival'.

- cex:

  Dimension of the points used for plotting the estimates. Defaults to
  0.2.

- cexlegend:

  Dimension of the text used for the legend. Defaults to 0.9.

## Value

The plot of the conditional survival function.

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

plot_survivalAdPaik(result)

 # } 
```
