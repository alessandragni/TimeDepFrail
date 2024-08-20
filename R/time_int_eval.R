#' @title Resolution of integral with respect to time
#'
#' @description
#' Function for the resolution of an integral with respect to time,
#' in the evaluation of the log-likelihood function.
#' It is implemented as defined in the paper of Wintrebert et al.'s (2004)
#'
#' @param time_t Event time instant
#' @param k k-th interval of the time-axis
#' @param time_axis Temporal domain (it may coincide with the follow-up)
#'
#' @return Evaluation of the temporal integral

time_int_eval <- function(time_t, k, time_axis){
  if (time_t < time_axis[k])
    return (0)
  else if (time_t >= time_axis[k] & time_t < time_axis[k+1])
    return (time_t - time_axis[k])
  else if(time_t >= time_axis[k+1])
    return (time_axis[k+1] - time_axis[k])
}
