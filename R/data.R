#' Data Dropout Dataset
#'
#' This dataset is extracted from an administrative database provided by a 
#' non-specified university and tracks students enrolled in 2012 over three academic years (or 6 semesters).
#'  We are interested in understanding what factors lead to students dropping out. 
#'  Dropout students with a time-instant in the first semester have been removed, 
#'  for internal reasons (the university cannot take preventive action to reduce or avoid their withdrawal).

#' The students are followed for at most 3 academic years or, equivalently, 
#' 6 semesters (follow-up periods), from the first day of lecture up to the time-instant 
#' of withdrawal (i.e. survival event) or the end of the academic year.
#'
#' @format A data frame with 4448 rows and 4 columns:
#' \describe{
#'   \item{Gender}{Categorical covariate (Male or Female).}
#'   \item{CFUP}{Standardized numerical covariate indicating the number of credits \ passed by the students by the end of the first semester.}
#'   \item{time_to_event}{Time (in semesters) at which a student decides to leave the university. \ A value greater than 6.0 indicates the student did not drop out during the follow-up (e.g. 6.1 semesters)}
#'   \item{group}{Categorical variable indicating the student's course of study, with 16 different levels from CosA, CosB, ... , CosP.}
#' }
#' @source Data for demonstration purposes.
"data_dropout"
