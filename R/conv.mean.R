#'
#' @title Running mean function
#'
#' @description A running mean, intended to check when convergence is achieved (useful for setting mcsize or bssize).
#'
#' @param x a variable whose convergence we want to check
#'
#' @return returns a running mean of the variable whose convergence we want to check
#' @export
#'
#' @import stats
#'
#'
conv.mean <- function(x) {cumsum(x)/(1:length(x))}
