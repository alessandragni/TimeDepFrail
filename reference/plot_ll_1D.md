# Plot the One-Dimensional Log-Likelihood Function

This function plots the trend of the log-likelihood function concerning
a single parameter specified by its index in the parameter vector. It
generates samples of the parameter, evaluates them in the log-likelihood
function, and displays the results along with the maximum point of the
one-dimensional log-likelihood function.

## Usage

``` r
plot_ll_1D(
  param_1D,
  index_param_1D,
  ll_1D,
  params,
  param_range_min,
  param_range_max,
  dataset,
  centre,
  time_axis,
  dropout_matrix,
  e_matrix,
  n_points = 150,
  cex = 0.7,
  cex_max = 0.8,
  color_bg = "black",
  color_max_bg = "red",
  pch = 21
)
```

## Arguments

- param_1D:

  A numeric value representing the optimal parameter determined by
  maximizing the log-likelihood function for the specified parameter.

- index_param_1D:

  An integer representing the index of the optimal parameter within the
  parameter vector.

- ll_1D:

  A numeric value of the log-likelihood function evaluated at the
  optimal parameter `param_1D`, with the other parameters held constant.

- params:

  A numeric vector of length equal to the number of parameters minus
  one, containing the fixed values for the other parameters.

- param_range_min:

  A numeric value indicating the minimum allowable value for the
  parameter `param_1D`.

- param_range_max:

  A numeric value indicating the maximum allowable value for the
  parameter `param_1D`.

- dataset:

  A data frame or matrix containing individual covariates.

- centre:

  A numeric vector indicating individual cluster membership; its length
  must match the number of individuals in the dataset.

- time_axis:

  A numeric vector corresponding to the subdivisions of the temporal
  domain.

- dropout_matrix:

  A binary matrix indicating which interval of the time domain an
  individual failed. Each row should sum to 1 (if failed) or 0 (if not
  failed), with dimensions (n_individuals, n_intervals).

- e_matrix:

  A matrix of dimensions (n_individuals, n_intervals), where each
  element contains the evaluation of the temporal integral performed by
  the function `time_int_eval`.

- n_points:

  An integer specifying the number of points at which to evaluate the
  log-likelihood function. A value that is neither too small nor too
  high is recommended; the default is 150.

- cex:

  A numeric value specifying the size of the points used for the
  graphical representation of the log-likelihood function. Default is
  0.7.

- cex_max:

  A numeric value indicating the size of the optimal point (the one
  maximizing the log-likelihood function). Default is 0.8.

- color_bg:

  A string specifying the color for the points representing the
  log-likelihood trend. Default is `'black'`.

- color_max_bg:

  A string specifying the color for the optimal point provided as the
  first argument. Default is `'red'`.

- pch:

  A numeric or character symbol representing the shape of the plotted
  points. Default is a circle (`21`).

## Value

A plot displaying the trend of the log-likelihood function concerning a
single parameter, including the maximum point.
