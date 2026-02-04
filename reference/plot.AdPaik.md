# Plots Related to the the 'Adapted Paik et Al.' Model

Plots Related to the the 'Adapted Paik et Al.' Model

## Usage

``` r
# S3 method for class 'AdPaik'
plot(
  x,
  which = c(1, 2),
  captions = c("Plot 1: Baseline Hazard", "Plot 2: Posterior Frailty Estimate"),
  ...
)
```

## Arguments

- x:

  An object of class 'AdPaik'.

- which:

  A numeric vector indicating which plots to display. Choices: 1 =
  Baseline Hazard, 2 = Posterior Frailty Estimate.

- captions:

  A character vector with captions for each plot.

- ...:

  Additional arguments to be passed to other methods.

## Value

No return value. This function generates plots.

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
#> Error in ll_AdPaik_centre_1D(x, index, params, dataset_centre, dropout_matrix_centre,     e_matrix_centre): dims [product 330] do not match the length of object [3300]

plot(result)
#> Error: object 'result' not found
# }  
```
