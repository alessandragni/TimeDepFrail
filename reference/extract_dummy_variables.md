# Transform Categorical Covariate into Dummy Variables

This function produces for a categorical variable of the dataset
(covariate) the associated dummy variables: for n levels of the
covariate, (n-1) dummy binary variables are generated. The chosen
reference value is the first one of the list of extracted levels and
cannot be changed (the first one in alphabetical order). Therefore, if
an individual has null value for all dummy variables, then his/her
belonging level is the reference one.

Each dummy variable has a name, corresponding to the name of the
covariate + name of the level.

## Usage

``` r
extract_dummy_variables(covariate, covariate_name)
```

## Arguments

- covariate:

  Categorical dataset covariate, with at least 2 levels.

- covariate_name:

  Name of the covariate, for assigning each dummy variable a proper
  name.

## Value

S3 object of class 'DummyData', composed of three elements. See details.

## Details

The S3 class object 'DummyData' contains the variables related to the
transformation of a single categorical covariate present in the dataset
into (n-1) binary covariates, stored in a matrix. To be precise:

- DummyMatrix: binary matrix of dimension (n_individuals, n-1), where
  each column corresponds to one level of the original categorical
  covariate.

- DummyVariablesName: categorical vector of length (n-1), reporting the
  names of the dummy variables and, therefore, the new name of each
  regressor.

- DummyVariablesNumber: number of dummy variables (n-1).
