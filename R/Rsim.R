#' @title Reference-free simplex R-learner
#' @description Ancillary function to obtain \code{Rsim} results
#' @param x Data matrix
#' @param y.adj \eqn{Y - m(X)}
#' @param Trt Treatment assignment indicator
#' @param x.test Test dataset
#' @param pi.hat Estimated propensity score for x
#' @param pi.test Estimated propensity score for x.test
#' @param W Simplex coordinates
#' @return A matrix of treatment effect mainly for optimal treatment recommendation
#'

Rsim <- function(x, y.adj, Trt, x.test, pi.hat, pi.test, W) {
  # transform data X into required shape:
  nobs = nrow(x)
  pobs = ncol(x)
  k = ncol(W)
  x.whole = cbind(1, x)
  x.new = sapply(seq(nobs), function(i){
    as.vector(outer(W[,Trt[i]] ,x.whole[i,]))/pi.hat[i,Trt[i]]
  })
  x.new = t(x.new)
  penalty_f = c(rep(0,k-1), rep(1, pobs*(k-1)))
  fit.tau = cv.glmnet(x.new, y.adj, family = "gaussian", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)

  x.test.whole = cbind(1, x.test)
  best.beta = coef(fit.tau,s="lambda.min")
  best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
  est.tau = x.test.whole %*% best.beta %*% W / pi.test
  return(est.tau)
}



