# Check Structure for the Parameters Confidence Interval

The function controls that the structure of the Parameters Confidence
Intervals coincides with the theoretical one.

## Usage

``` r
check.structure_paramsCI(parametersCI)
```

## Arguments

- parametersCI:

  S3 object of class 'ParametersCI', composed of two elements:

  - left confidence interval: numerical vector of length equal to the
    number of parameters in the model

  - right confidence interval: numerical vector of length equal to the
    number of parameters in the model

## Value

An error if any condition is not satisfied.
