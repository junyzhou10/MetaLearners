#' @title Meta-learners for Treatment Recommendation and Individual Treatment Effect Estimation for multiple-treatment scenario
#' @description ITE/HTE/CATE estimation can be achieved by S-, T-, X-, R-, deC-learning,
#' while optimal treatment can be achieved by S-, T-, Rsim-, deC- and AD-learning.
#' @param X Data matrix
#' @param Y Outcome
#' @param Trt Treatment assignment indicator
#' @param X.test Test data. If \code{NULL}, training dataset, i.e., \code{X} will be used. In practice, \code{X.test} is always \code{NULL}.
#' @param Learners The S-, T-, and deC-learners are always adopted since they can do both treatment effect
#' estimation and optimal treatment recommendation. Other learners are optional. Default includes all of them.
#' @param algorithm Machine learning algorithm for the analysis. Currently, it supports one of \code{c("BART", "GAM", "RF", "SL")}.
#'               representing Bayesian Additive Regression Trees (BART), Generalized Boosted Method (GBM), Random Forest (RF), and Super Learner (SL).
#' @param controls Additional arguments for plot
#' @examples
#'
#'
#' @export

MetaLearners <- function(X,
                         Y,
                         Trt,
                         X.test = NULL,
                         Learners = c("X", "R", "Rsim", "AD"),
                         algorithm = "BART",
                         controls = list(SL.library = c("SL.bartMachine", "SL.gam", "SL.ranger"))
) {
  print("h")
}


