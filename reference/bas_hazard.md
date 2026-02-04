# Baseline Hazard Step-Function

The method computes the baseline hazard step-function in each interval
of the time-domain, using the estimated parameters \\\phi_k, \forall k\\

## Usage

``` r
bas_hazard(object)
```

## Arguments

- object:

  S3 object of class 'AdPaik' returned by the main model output, that
  contains all the information for the computation of the frailty
  standard deviation.

## Value

Numerical vector of length equal to the number of intervals of the
time-domain, with the value of the baseline hazard step-function.

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

bas_hazard(result)
#> Error: object 'result' not found
# }
```
