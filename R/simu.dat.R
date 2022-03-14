#' One simulated data with four treatment options
#'
#' A total of 400 subject receiving 4 treatments randomly but with unequal chance, the allocation is 1:1:1:7
#'
#' \itemize{
#'   \item X Observed covariates
#'   \item Y Observed continuous outcomes
#'   \item Trt Observed treatment indicators
#'   \item tau True treatment effect which can be used to evaluate the treatment effect and optimal treatment
#' }
#' @docType data
#' @name simu.dat
#' @format A list of 4 items. X: 400x30 matrix indicating 400 subjects and 30 covariate (but only part of them is related)
#'
"simu.dat"
