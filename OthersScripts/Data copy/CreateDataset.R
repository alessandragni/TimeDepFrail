# WINDOWS SETUP
# Set working directory, save the functions and load the dataset
setwd("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R")
load("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R/Data/data_dropout.RData")
#load("C:/Users/admin/Documents/R/Thesis/TimeDependentSharedFrailtyCoxModels-R/Data/dataless_time_varying_year2012.RData")

# MAC SETUP
setwd("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R")
#load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/data_dropout.RData")
#save.image("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/data_dropout.RData")
load("/Users/noa/Desktop/Giulia/TimeDependentSharedFrailtyCoxModels-R/Data/dataless_time_varying_year2012.RData")

setwd("~/Documents/DATA/POLITECNICO/PHD/CODE_REPO/TimeDepFrail/Data")
load("~/Documents/DATA/POLITECNICO/PHD/CODE_REPO/TimeDepFrail/Data/dataless_time_varying_year2012.RData")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# This script is only used to fix and change the dataset variables
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Remove useless variables
rm(cox_data_year, new_data_year, nodes_ghqm, nodesG_ghqm, weights_ghqm, weightsG_ghqm, a_interval)
rm(index_year, L, n_interval, n_node, n_nodeG, n_other_variable, n_regressor, year3)

# Change name for the course of study
old_name <- c("EngA", "EngB", "EngC", "EngD", "EngE","EngF", "EngG", "EngH", "EngI", "EngJ", "EngK", "EngL",
              "EngM", "EngN", "EngO", "EngP")
new_name <- c("CosA", "CosB", "CosC", "CosD", "CosE", "CosF", "CosG", "CosH", "CosI", "CosJ", "CosK", "CosL",
              "CosM", "CosN", "CosO", "CosP")
for(j in 1:length(old_name)){
  index <- which(faculty_codes == old_name[j])
  faculty_codes[index] <- new_name[j]
}

# Change levels for Gender
old_levels <- c(0,1)
new_levels <- c("Male", "Female")
for(j in 1:length(old_levels)){
  index <- which(data_app_year[,1] == old_levels[j])
  data_app_year[index,1] <- new_levels[j]
}

# Create new dataset
data_dropout <- data.frame(data_app_year)
data_dropout <- cbind(data_dropout, faculty_codes)
colnames(data_dropout) <- c("Gender", "CFUP", "time_to_event", "group")

# Round elements
second_column <- as.numeric(data_dropout[,2])
round(second_column, 7)
data_dropout[,2] <- second_column

# Remove useless variables
rm(faculty_codes, data_app_year, old_name, new_name)
rm(j, index, second_column)
rm(old_levels, new_levels)
# rm(element,  n_individuals)

#-------------------------------------------------------------------------------
save(data_dropout, file = "data_dropout.RData")

# # Save dataset into R package
# library(devtools)
# use_data(data_dropout, overwrite= TRUE)
#
# # Add documentation for data
# use_r("data")

