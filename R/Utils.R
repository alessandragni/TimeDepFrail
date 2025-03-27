format_perc <- function(probs, digits) {
  paste0(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), " %")
}


#' @title Resolution of Integral with Respect to Time
#'
#' @description
#' Function for the resolution of an integral with respect to time,
#' in the evaluation of the log-likelihood function.
#' It is implemented as defined in the paper of Wintrebert et al.'s (2004)
#'
#' @param time_t Event time instant
#' @param k k-th interval of the time-axis
#' @param time_axis Temporal domain
#'
#' @return Evaluation of the temporal integral
#' 
#' @keywords internal

time_int_eval <- function(time_t, k, time_axis){
  if (time_t < time_axis[k])
    return (0)
  else if (time_t >= time_axis[k] & time_t < time_axis[k+1])
    return (time_t - time_axis[k])
  else if(time_t >= time_axis[k+1])
    return (time_axis[k+1] - time_axis[k])
}
