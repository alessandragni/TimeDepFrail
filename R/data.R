#' Data Dropout Dataset
#'
#' This dataset contains information about dropout rates in a study.
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
