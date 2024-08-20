#' @title
#' Transform categorical covariate into dummy variables
#'
#' @description
#' This function produces for a categorical variable of the dataset (covariate) the associated dummy variables: for n levels of the covariate,
#' (n-1) dummy binary variables are generated. The chosen reference value is the first one of the list of extracted levels and cannot be changed
#' (the first one in alphabetical order). Therefore, if an individual has null value for all dummy variables, then his/her belonging level is the
#' reference one.
#'
#' Each dummy variable has a name, corresponding to the name of the covariate + name of the level.
#'
#'
#' @param covariate Categorical dataset covariate, with at least 2 levels.
#' @param covariate_name Name of the covariate, for assigning each dummy variable a proper name.
#'
#' @return S3 object of class 'DummyData', composed of three elements. See details.
#'
#' @details
#' The S3 class object 'DummyData' contains the variables related to the transformation of a single categorical covariate present in the dataset
#' into (n-1) binary covariates, stored in a matrix. To be precise:
#' - DummyMatrix: binary matrix of dimension (n_individuals, n-1), where each column corresponds to one level of the original categorical
#' covariate.
#' - DummyVariablesName: categorical vector of length (n-1), reporting the names of the dummy variables and, therefore, the new name of each regressor.
#' - DummyVariablesNumber: number of dummy variables (n-1).

extract_dummy_variables <- function(covariate, covariate_name){
  # Extract number of individuals
  n_individuals <- length(covariate)
  
  #Extract levels from covariate vector
  levels_covariate <- levels(factor(covariate))
  
  # Compute the number of levels
  n_levels_covariate <- length(levels_covariate)
  
  # Check correctness of covariate levels numerosity
  if(n_levels_covariate == 0){
    error_msg <- paste("No levels for covariate", covariate_name)
    stop(error_msg)
  }
  else if(n_levels_covariate == 1){
    error_msg <- paste("Only one level for covariate", covariate_name)
    stop(error_msg)
  }
  
  # Create matrix for dummy variables
  n_dummy_variables <- n_levels_covariate - 1
  dummy_matrix <- matrix(rep(0, n_individuals * n_dummy_variables), n_individuals, n_dummy_variables)
  dummy_name <- rep("", n_dummy_variables)
  
  # Select individuals belonging to levels (d+1), so that reference one is (1)
  for(d in 1:n_dummy_variables){
    # Extract individuals in level (d+1)
    individuals_level_d <- which(covariate == levels_covariate[d+1])
    dummy_matrix[individuals_level_d, d] <- 1
    
    # Create also name for the matrix columns
    column_name <- paste(covariate_name, levels_covariate[d+1], sep="", collapse="")
    dummy_name[d] <- column_name
  }
  colnames(dummy_matrix) <- dummy_name
  
  # Return dummy matrix and columns name
  return_list <- list("DummyMatrix" = dummy_matrix,
                      "DummyVariablesName" = dummy_name,
                      "DummyVariablesNumber" = n_dummy_variables)
  class(return_list) <- "DummyData"
  
  return(return_list)
}
