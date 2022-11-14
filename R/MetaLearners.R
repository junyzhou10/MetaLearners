#' @title Meta-learners for Treatment Recommendation and Individual Treatment Effect Estimation for multiple-treatment scenario
#' @description ITE/HTE/CATE estimation can be achieved by S-, T-, X-, R-, deC-learning,
#' while optimal treatment can be achieved by S-, T-, Rsim-, deC- and AD-learning.
#' @param X Data matrix
#' @param Y Outcome. Support either continuous or binary. If \code{Y} only has two levels, will be deemed
#'          as binary outcome
#' @param Trt Treatment assignment indicator
#' @param X.test Test data. If \code{NULL}, training dataset, i.e., \code{X} will be used. In practice, \code{X.test} is always \code{NULL}.
#' @param Learners The S-, T-, and deC-learners are always adopted since they can do both treatment effect
#' estimation and optimal treatment recommendation. Other learners are optional. Default includes all of them.
#' @param algorithm Machine learning algorithm for the analysis. Currently, it supports one of \code{c("BART", "GAM", "RF", "SL")}.
#'               representing Bayesian Additive Regression Trees (BART), Generalized Boosted Method (GBM), Random Forest (RF),
#'               and Super Learner (SL). The default is \code{BART}.
#' @param controls A list of control arguments for analysis. See details.
#'
#' @details For S-, T-, X-, and R-learner, they are initially proposed for two-treatment scenarios,
#' but can be extended to multiple-treatment following the same gist. R-learner will yield multiple sets
#' of results with the choise of reference group, and the results can be different from each other,
#' especially when treatment allocation is unbalanced. This algorithm will return results from all possible
#' reference group, named as \code{R.k}. User can then decide which reference group to take (technically, the group with largest
#' sample is more suitable as reference one). \cr\cr
#' R-, simplex R (Rsim), deC-, and AD-learning have to choose the basis function since a global additive structure
#' is assumed. Different basis functions are available and user can adjust the corresponding \code{controls} arguments.
#' Interaction terms are not available to generate automatically but user can input their own design matrix.
#' The lasso-type penalization from \code{glmnet} will be adopted in analysis for R-, simplex R (Rsim), and deC-learner,
#' and group lasso for AD-learning.
#' \cr\cr
#' \code{controls} is a list containing:
#' \itemize{
#' \item{SL.library}{: The algorithms to be included in for Super Learner. For algorithm is \code{SL} only.
#' For details, please check package [SuperLearner]}
#' \item{SL.library.PS}{: The SL.library for propensity score estimation. It can be different from \code{SL.library}
#' given some algorithms may not capable for multinomial outcomes. For algorithm is \code{SL} only.}
#' \item{basis.func}{: To specify the basis functions. Currently support one of these
#' \code{c("Polynomial", "bs", "ns")}, standing for polynomial regression, B-splines, and natural cubic spline.}
#' \item{degree}{: The degree for basis functions (not for natural cubic spline)}
#' \item{n.knots}{: Number of knots when B-slines or natural spline is adopted.}
#' \item{X_train}{: User provided design matrix. If given, it overrides all basis function related inputs.}
#' \item{X_test}{: Similar to \code{X_train}. For really world application rather than simulation studies,
#' it is not needed.}
#' }
#'
#' @return
#' \item{S.res}{Results from S-learner. If binary outcome, then it represents the probability of \code{T=1}.
#'              Causal estimands like log(Relative Risk), or log(Odds Ratio) can be obtained
#'              from these estimated probabilities. Same for the rest.}
#' \item{T.res}{Results from T-learner}
#' \item{X.res}{Results from X-learner}
#' \item{R.res}{Results from R-learner. Notably, it contains multiple list of results with different
#' choice of reference. The results can be different with the reference level selection, especially
#' when treatment allocation is unbalanced.}
#' \item{Rsim.res}{Results from simplex R-learner}
#' \item{C.res}{Results from deC-learner, which contains \code{C.resS}, \code{C.resT}, and \code{C.resST},
#' representing using S-, T-, and average of S- and T- in calculating the model average. When sample size
#' is small or treatment allocation is unbalanced, \code{C.resS} is more preferred.}
#' \item{AD.res}{Results from AD-learning}
#' @examples
#' \dontrun{
#' require(doMC)
#' registerDoMC(cores = 6)
#' res = MetaLearners(simu.dat$X, simu.dat$Y, simu.dat$Trt, algorithm = "GAM")
#' final.res = parseRes(res)
#' }
#' @md
#' @import SuperLearner ranger splines dbarts BART glmnet dplyr
#' @importFrom stats coef predict pnorm
#' @export

MetaLearners <- function(X,
                         Y,
                         Trt,
                         X.test = NULL,
                         Learners = c("X", "R", "Rsim", "AD"),
                         algorithm = "BART",
                         controls = list(
                           SL.library = c("SL.bartMachine", "SL.gam", "SL.ranger"),
                           SL.library.PS = c("SL.bartMachine", "SL.gam", "SL.ranger"),
                           basis.func = "Polynomial",
                           degree = 2,
                           n.knots = 3,
                           X_train = NULL,
                           X_test = NULL
                         )
) {
  # test dataset
  if (is.null(X.test)) {
    X.test = X
  }

  # detect number of treatments & generate global variables
  K.grp = sort(unique(Trt))
  k = length(K.grp)
  p = ncol(X)
  n.train = nrow(X)
  n.test  = nrow(X.test)

  if (k<=2) {
    stop("The algorithm is designed for multiple treatment, currently only support for k>2.")
  }
  # check algorithm
  if (length(algorithm)!=1 | !algorithm %in% c("BART", "GAM", "RF", "SL") ){
    stop("The algorithm at this moment only support one of the following inputs for `algorithm`: 'BART', 'GAM', 'RF', 'SL'!")
  }


  # check controls:
  controls.default = list(
    SL.library = c("SL.bartMachine", "SL.gam", "SL.ranger"),
    SL.library.PS = c("SL.bartMachine", "SL.gam", "SL.ranger"),
    basis.func = "Polynomial",
    degree = 2,
    n.knots = 3,
    X_train = NULL,
    X_test = NULL
  )
  name.ctrl = names(controls.default)
  for (ii in seq(length(name.ctrl))) {
    if (is.null(controls[[name.ctrl[ii]]])) {
      controls[[name.ctrl[ii]]] <- controls.default[[name.ctrl[ii]]]
    }
  }

  # generate simplex coordinates
  k = length(K.grp)
  z = rep(1, k-1)
  e = diag(x = 1, k-1)
  W = cbind((k-1)^(-0.5) * z,  (k/(k-1))^(0.5)*e - z*(1+sqrt(k))/(k-1)^1.5)

  # generate design matrix following arguments in `controls`
  make_matrix = function(x) stats::model.matrix(~.-1, x)
  if (is.null(controls$X_train)) {
    ### Polynomial
    if (controls$basis.func == "Polynomial") {
      X_bs = NULL; x.tmp = rbind(X, X.test)
      for (nn in seq(controls$degree)) {
        X_bs = cbind(X_bs, x.tmp^nn)
      }
      X_bs = apply(X_bs, 2, function(x) scale(x, scale = F)) # avoid collinearity
    }
    ### B-spline
    if (controls$basis.func == "bs") {
      X_bs = do.call(cbind, lapply(1:p, function(col){matrix(
        splines::bs(rbind(X, X.test)[, col],
                    df = controls$n.knots + controls$degree,
                    degree = controls$degree, intercept = F),
        n.train+n.test)}))
    }
    ### natural cubic spline
    if (controls$basis.func == "ns") {
      X_bs = do.call(cbind, lapply(1:p, function(col){matrix(
        splines::ns(rbind(X, X.test)[, col],
                    df = controls$n.knots + 2,
                    intercept = F),
        n.train+n.test)}))
    }
    # make design matrix
    X_train = data.frame(X_bs[1:n.train,]) %>% make_matrix
    X_test = data.frame(X_bs[-c(1:n.train),]) %>% make_matrix
  } else { # design matrix is provided
    X_train = controls$X_train
    if (is.null(controls$X_test)) {
      X_test = X_train
    } else {
      X_test = controls$X_test
    }
    X_train = data.frame(X_train) %>% make_matrix
    X_test = data.frame(X_test) %>% make_matrix
  }

  pobs = ncol(X_train)
  # returns
  S.res <- T.res <- X.res <- R.res <- Rsim.res <- C.res <- AD.res <- NULL


  ############ Diverge by outcome type:
  if (length(unique(Y)) == 2) { # Binary
    ######=============================== BART ==================================#####
    if (algorithm == "BART") {
      ######  BART: S-learner  ######
      dat.train = data.frame(trt = Trt, X)
      dat.tmp = dat.train
      for (ii in K.grp) {
        tmp = data.frame(trt = ii, rbind(X, X.test))
        dat.tmp = rbind(dat.tmp, tmp)
      }
      dat.tmp[,1] = as.factor(dat.tmp[,1])
      Z = as.matrix(dat.tmp[,-1])
      dat.S = stats::model.matrix(~trt*Z-1, dat.tmp)
      dat.train = dat.S[1:n.train,]; dat.test = dat.S[-(1:n.train),]
      S.fit = bart(x.train = dat.train, y.train = Y, x.test = dat.test, ntree = 200, verbose = FALSE)
      S.res = matrix(colMeans(S.fit$yhat.test), ncol = length(K.grp), byrow = F)
      S.res = pnorm(S.res) # Bart use probit link
      SC.res= S.res[1:n.train,]; SC.res = log(SC.res/(1-SC.res))  # for deC learner, we need logOR, so here we report log(Odds)
      S.res = S.res[-(1:n.train),] # here we report is mu(x) = Pr(Y=1; x)


      ######  BART: T-learner  ######
      T.res = NULL; TX.res = NULL
      for (ii in K.grp) {
        T.bart = bart(x.train = X[Trt==ii,], y.train = Y[Trt==ii], x.test = rbind(X, X.test), ntree = 200, verbose = FALSE)
        TX.res = cbind(TX.res, colMeans(T.bart$yhat.test)[1:nrow(X)])
        T.res = cbind(T.res, colMeans(T.bart$yhat.test)[-(1:nrow(X))])
      }
      T.res = pnorm(T.res) # here we report is mu(x) = Pr(Y=1; x)
      TX.res = pnorm(TX.res); TX.res = log(TX.res/(1-TX.res)) # log(Odds) for deC-learner

      ######  BART: X-learner  ######
      # not available so far

      ######  BART: R/Rsim-learner  ######
      # not available for Robinson's decomposition
    }


    ######=============================== GAM ==================================#####
    if (algorithm == "GAM") {
      ######  GAM: S-learner  ######
      dat.train = data.frame(trt = Trt, X_train)
      dat.tmp = dat.train
      for (ii in K.grp) {
        tmp = data.frame(trt = ii, rbind(X_train, X_test))
        dat.tmp = rbind(dat.tmp, tmp)
      }
      dat.tmp[,1] = as.factor(dat.tmp[,1])
      Z = as.matrix(dat.tmp[,-1])
      dat.S = stats::model.matrix(~trt*Z-1, dat.tmp)
      dat.train = dat.S[1:n.train,]; dat.test = dat.S[-(1:n.train),]
      S.GAM = cv.glmnet(dat.train, Y, family = "binomial", parallel = TRUE, maxit = 100000, intercept=T)
      S.res = matrix(predict(S.GAM, newx = dat.test, type = "response", s = "lambda.min"), ncol = length(K.grp), byrow = F)
      SC.res= S.res[1:n.train,]; SC.res = log(SC.res/(1-SC.res))
      S.res = S.res[-(1:n.train),]

      ######  GAM: T-learner  ######
      T.res = NULL; TX.res = NULL
      for (ii in K.grp) {
        T.GAM = cv.glmnet(X_train[Trt==ii,], Y[Trt==ii], family = "binomial", parallel = TRUE, maxit = 100000, intercept=T)
        TX.res = cbind(TX.res, predict(T.GAM, newx = X_train, type = "response", s = "lambda.min"))
        T.res = cbind(T.res, predict(T.GAM, newx = X_test, type = "response", s = "lambda.min"))
      }
      TX.res = log(TX.res/(1-TX.res)) # log(Odds) for deC-learner

    }


    if (algorithm == "RF") {
      ######  RF: S-learner  ######
      dat.train = data.frame(trt = Trt, X)
      dat.tmp = dat.train
      for (ii in K.grp) {
        tmp = data.frame(trt = ii, rbind(X, X.test))
        dat.tmp = rbind(dat.tmp, tmp)
      }
      dat.tmp[,1] = as.factor(dat.tmp[,1])
      Z = as.matrix(dat.tmp[,-1])
      dat.S = stats::model.matrix(~trt*Z-1, dat.tmp)
      dat.train = dat.S[1:n.train,]; dat.test = dat.S[-(1:n.train),]
      colnames(dat.test) <- colnames(data.frame(dat.train))
      S.RF = ranger(Y~., data = data.frame(Y = Y, dat.train), num.trees = 500, probability = T)
      S.res = matrix(predict(S.RF, data = dat.test)$predictions[,1], ncol = length(K.grp), byrow = F)
      SC.res= S.res[1:n.train,]; SC.res = log(SC.res/(1-SC.res))
      S.res = S.res[-(1:n.train),]

      ######  RF: T-learner  ######
      T.res = NULL; TX.res = NULL
      for (ii in K.grp) {
        T.RF = ranger(Y~., data = data.frame(Y = Y[Trt==ii], X[Trt==ii,]), num.trees = 500, probability = T)
        TX.res = cbind(TX.res, predict(T.RF, data = data.frame(X))$predictions[,1])
        T.res = cbind(T.res, predict(T.RF, data = data.frame(X.test))$predictions[,1])
      }
      TX.res = log(TX.res/(1-TX.res)) # log(Odds) for deC-learner
    }


    ######=================  General: deC-learner  =====================######
    h.hat  = (rowMeans(TX.res) + rowMeans(SC.res))/2
    h.hatS = rowMeans(SC.res)
    h.hatT = rowMeans(TX.res)

    # transform data X into required shape:
    x.whole = cbind(1, X_train)
    x.new = sapply(seq(n.train), function(i){
      as.vector(outer(W[,Trt[i]] ,x.whole[i,]))
    })
    x.new = t(x.new)
    penalty_f = c(rep(0,k-1), rep(1, pobs*(k-1)))
    fit.tau  = cv.glmnet(x = x.new, y = Y, offset = h.hat, family = "binomial", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)
    fit.tauS = cv.glmnet(x = x.new, y = Y, offset = h.hatS, family = "binomial", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)
    fit.tauT = cv.glmnet(x = x.new, y = Y, offset = h.hatT, family = "binomial", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)

    x.test.whole = cbind(1, X_test)
    # combine both S- and T- results
    best.beta = stats::coef(fit.tau,s="lambda.min")
    best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
    C.resST = x.test.whole %*% best.beta %*% W + h.hat # this is log(Odds)
    # based on S-
    best.beta = stats::coef(fit.tauS, s="lambda.min")
    best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
    C.resS = x.test.whole %*% best.beta %*% W + h.hatS
    # based on T-
    best.beta = stats::coef(fit.tauT, s="lambda.min")
    best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
    C.resT = x.test.whole %*% best.beta %*% W + h.hatT
    # deC-learner results
    C.res = list(C.resST = 1/(1+exp(-C.resST)), C.resS = 1/(1+exp(-C.resS)), C.resT = 1/(1+exp(-C.resT))) # keep consistent, return Pr(Y=1)

    return(list(S.res = S.res, T.res = T.res, C.res = C.res))

  } else { # continuous


    ######=============================== BART ==================================#####
    if (algorithm == "BART") {
      ######  BART: S-learner  ######
      dat.train = data.frame(trt = Trt, X)
      dat.tmp = dat.train
      for (ii in K.grp) {
        tmp = data.frame(trt = ii, rbind(X, X.test))
        dat.tmp = rbind(dat.tmp, tmp)
      }
      dat.tmp[,1] = as.factor(dat.tmp[,1])
      Z = as.matrix(dat.tmp[,-1])
      dat.S = stats::model.matrix(~trt*Z-1, dat.tmp)
      dat.train = dat.S[1:n.train,]; dat.test = dat.S[-(1:n.train),]
      S.fit = bart(x.train = dat.train, y.train = Y, x.test = dat.test, ntree = 200, verbose = FALSE)
      S.res = matrix(colMeans(S.fit$yhat.test), ncol = length(K.grp), byrow = F)
      SC.res= S.res[1:n.train,]; S.res = S.res[-(1:n.train),]


      ######  BART: T-learner  ######
      T.res = NULL; TX.res = NULL
      for (ii in K.grp) {
        T.bart = bart(x.train = X[Trt==ii,], y.train = Y[Trt==ii], x.test = rbind(X, X.test), ntree = 200, verbose = FALSE)
        TX.res = cbind(TX.res, T.bart$yhat.test.mean[1:nrow(X)])
        T.res = cbind(T.res, T.bart$yhat.test.mean[-(1:nrow(X))])
      }

      # check if propensity score estimation is required
      if (sum(Learners %in% c("X","R","Rsim","AD")) > 0) {
        invisible(utils::capture.output(PS.bart <- mbart(x.train = X, y.train = Trt, x.test = rbind(X,X.test), ntree = 200)))
        ps.bart = matrix(PS.bart$prob.test.mean, ncol = length(K.grp), byrow = T)
        ps.train = ps.bart[1:nrow(X),]
        ps.test  = ps.bart[-(1:nrow(X)),]
      }

      ######  BART: X-learner  ######
      if ("X" %in% Learners) {
        X.res = NULL; X.names = NULL
        # Imputed outcome responses (all potential outcomes are imputed):
        for (i in K.grp) { # Y_i - \hat{Y}_j & \hat{Y}_i - Y_j
          for (j in K.grp[K.grp > i]) {
            X.bart.i = bart(x.train = X[Trt==i,], y.train = Y[Trt==i] - TX.res[Trt==i,j], x.test = X.test, ntree = 200, verbose = FALSE)
            X.bart.j = bart(x.train = X[Trt==j,], y.train = TX.res[Trt==j,i] - Y[Trt==j], x.test = X.test, ntree = 200, verbose = FALSE)
            ps.ij    = ps.test[,c(i,j)]
            ps.ij    = ps.ij/rowSums(ps.ij)
            X.res    = cbind(X.res, X.bart.i$yhat.test.mean*ps.ij[,2]+X.bart.j$yhat.test.mean*ps.ij[,1])
            X.names  = c(X.names, paste0(i,"-",j))
          }
        }
        colnames(X.res) <- X.names
      }


      ######  BART: R/Rsim-learner  ######
      if ("R" %in% Learners | "Rsim" %in% Learners) {
        # estimate m(X) := E(Y|X)
        m.fit = bart(x.train = X, y.train = Y, x.test = X, ntree = 200, verbose = FALSE)
        est.m = m.fit$yhat.test.mean
        y.adj = Y - est.m
      }
    }

    ######=============================== GAM ==================================#####
    if (algorithm == "GAM") {
      ######  GAM: S-learner  ######
      dat.train = data.frame(trt = Trt, X_train)
      dat.tmp = dat.train
      for (ii in K.grp) {
        tmp = data.frame(trt = ii, rbind(X_train, X_test))
        dat.tmp = rbind(dat.tmp, tmp)
      }
      dat.tmp[,1] = as.factor(dat.tmp[,1])
      Z = as.matrix(dat.tmp[,-1])
      dat.S = stats::model.matrix(~trt*Z-1, dat.tmp)
      dat.train = dat.S[1:n.train,]; dat.test = dat.S[-(1:n.train),]
      S.GAM = cv.glmnet(dat.train, Y, family = "gaussian", parallel = TRUE, maxit = 100000, intercept=T)
      S.res = matrix(predict(S.GAM, newx = dat.test, type = "response", s = "lambda.min"), ncol = length(K.grp), byrow = F)
      SC.res= S.res[1:n.train,]; S.res = S.res[-(1:n.train),]

      ######  GAM: T-learner  ######
      T.res = NULL; TX.res = NULL
      for (ii in K.grp) {
        T.GAM = cv.glmnet(X_train[Trt==ii,], Y[Trt==ii], family = "gaussian", parallel = TRUE, maxit = 100000, intercept=T)
        TX.res = cbind(TX.res, predict(T.GAM, newx = X_train, type = "response", s = "lambda.min"))
        T.res = cbind(T.res, predict(T.GAM, newx = X_test, type = "response", s = "lambda.min"))
      }

      # check if propensity score estimation is required
      if (sum(Learners %in% c("X","R","Rsim","AD")) > 0) {
        if ( length(K.grp) > 2 ) {
          fit.pi  = cv.glmnet(X_train, Trt, family = "multinomial", type.measure = "deviance",parallel = TRUE, maxit = 100000, intercept=T)
          pi.hat  = predict(fit.pi, newx = X_test, s = "lambda.min", type = 'response')
          ps.test = pi.hat[,,1]
          pi.hat  = predict(fit.pi, newx = X_train, s = "lambda.min", type = 'response')
          ps.train= pi.hat[,,1]
        } else {
          fit.pi  = cv.glmnet(X_train, Trt, family = "binomial", type.measure = "deviance", parallel = TRUE, maxit = 100000, intercept=T)
          pi.hat  = predict(fit.pi, newx = X_test, s = "lambda.min", type = 'response')
          ps.test = cbind(1-pi.hat, pi.hat) + 1
          pi.hat  = predict(fit.pi, newx = X_train, s = "lambda.min", type = 'response')
          ps.train= cbind(1-pi.hat, pi.hat) + 1
          colnames(ps.test) <- colnames(ps.train) <- c("0", "1")
        }
      }

      ######  GAM: X-learner  ######
      if ("X" %in% Learners) {
        XX.res = NULL; X.names = NULL
        # Imputed outcome responses (all potential outcomes are imputed):
        for (i in K.grp) { # Y_i - \hat{Y}_j & \hat{Y}_i - Y_j
          for (j in K.grp[K.grp > i]) {
            X.gam.i = cv.glmnet(X_train[Trt==i,], Y[Trt==i] - TX.res[Trt==i,j], family = "gaussian", parallel = TRUE, maxit = 100000, intercept=T)
            X.gam.j = cv.glmnet(X_train[Trt==j,], TX.res[Trt==j,i] - Y[Trt==j], family = "gaussian", parallel = TRUE, maxit = 100000, intercept=T)
            ps.ij    = ps.test[,c(i,j)]
            ps.ij    = ps.ij/rowSums(ps.ij)
            X.res    = cbind(X.res, predict(X.gam.i, newx = X_test, type = "response", s = "lambda.min")*ps.ij[,2] + predict(X.gam.j, newx = X_test, type = "response", s = "lambda.min")*ps.ij[,1])
            X.names  = c(X.names, paste0(i,"-",j))
          }
        }
        colnames(X.res) <- X.names
      }


      ######  GAM: R/Rsim-learner  ######
      if ("R" %in% Learners | "Rsim" %in% Learners) {
        # estimate m(X) := E(Y|X)
        penalty_factor_nuisance = rep(1, pobs)
        fit.m = cv.glmnet(X_train, Y, family = "gaussian", parallel = TRUE, maxit = 100000, penalty.factor = penalty_factor_nuisance)
        est.m = predict(fit.m, newx = X_train, type = "response", s = "lambda.min")
        y.adj = Y - est.m
      }

    }

    ######=============================== RF ==================================#####
    if (algorithm == "RF") {
      ######  RF: S-learner  ######
      dat.train = data.frame(trt = Trt, X)
      dat.tmp = dat.train
      for (ii in K.grp) {
        tmp = data.frame(trt = ii, rbind(X, X.test))
        dat.tmp = rbind(dat.tmp, tmp)
      }
      dat.tmp[,1] = as.factor(dat.tmp[,1])
      Z = as.matrix(dat.tmp[,-1])
      dat.S = stats::model.matrix(~trt*Z-1, dat.tmp)
      dat.train = dat.S[1:n.train,]; dat.test = dat.S[-(1:n.train),]
      colnames(dat.test) <- colnames(data.frame(dat.train))
      S.RF = ranger(Y~., data = data.frame(Y = Y, dat.train), num.trees = 500)
      S.res = matrix(predict(S.RF, data = dat.test)$predictions, ncol = length(K.grp), byrow = F)
      SC.res= S.res[1:n.train,]; S.res = S.res[-(1:n.train),]

      ######  RF: T-learner  ######
      T.res = NULL; TX.res = NULL
      for (ii in K.grp) {
        T.RF = ranger(Y~., data = data.frame(Y = Y[Trt==ii], X[Trt==ii,]), num.trees = 500)
        TX.res = cbind(TX.res, predict(T.RF, data = data.frame(X))$predictions)
        T.res = cbind(T.res, predict(T.RF, data = data.frame(X.test))$predictions)
      }

      # check if propensity score estimation is required
      if (sum(Learners %in% c("X","R","Rsim","AD")) > 0) {
        # calculate propensity score
        PS.RF   = ranger(Trt ~., data = data.frame(Trt = as.factor(Trt), X), num.trees = 500, probability = T)
        ps.test = predict(PS.RF, data = data.frame(X.test))$predictions
        ps.train= predict(PS.RF, data = data.frame(X))$predictions
      }

      ######  RF: X-learner  ######
      if ("X" %in% Learners) {
        X.res = NULL; X.names = NULL
        # Imputed outcome responses (all potential outcomes are imputed):
        for (i in K.grp) { # Y_i - \hat{Y}_j & \hat{Y}_i - Y_j
          for (j in K.grp[K.grp > i]) {
            X.rf.i = ranger(tau ~., data = data.frame(tau = Y[Trt==i] - TX.res[Trt==i,j], X[Trt==i,]), num.trees = 500)
            X.rf.j = ranger(tau ~., data = data.frame(tau = TX.res[Trt==j,i] - Y[Trt==j], X[Trt==j,]), num.trees = 500)
            ps.ij    = ps.test[,c(i,j)]
            ps.ij    = ps.ij/rowSums(ps.ij)
            X.res    = cbind(X.res, predict(X.rf.i, data = data.frame(X.test))$predictions*ps.ij[,2]+predict(X.rf.j, data = data.frame(X.test))$predictions*ps.ij[,1])
            X.names  = c(X.names, paste0(i,"-",j))
          }
        }
        colnames(X.res) <- X.names
      }


      ######  RF: R/Rsim-learner  ######
      if ("R" %in% Learners | "Rsim" %in% Learners) {
        # estimate m(X) := E(Y|X)
        m.fit = ranger(Y~., data = data.frame(Y=Y, X), num.trees = 500)
        est.m = m.fit$predictions
        y.adj = Y - est.m
      }
    }

    ######=============================== SL ==================================#####
    if (algorithm == "SL") {
      stop("Still under development.")
    }


    ######=================  General: deC-learner  =====================######
    h.hat  = (rowMeans(TX.res) + rowMeans(SC.res))/2
    h.hatS = rowMeans(SC.res)
    h.hatT = rowMeans(TX.res)

    # transform data X into required shape:
    x.whole = cbind(1, X_train)
    x.new = sapply(seq(n.train), function(i){
      as.vector(outer(W[,Trt[i]] ,x.whole[i,]))
    })
    x.new = t(x.new)
    penalty_f = c(rep(0,k-1), rep(1, pobs*(k-1)))
    fit.tau  = cv.glmnet(x.new, Y-h.hat, family = "gaussian", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)
    fit.tauS = cv.glmnet(x.new, Y-h.hatS, family = "gaussian", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)
    fit.tauT = cv.glmnet(x.new, Y-h.hatT, family = "gaussian", parallel = TRUE, maxit = 100000, penalty.factor = penalty_f, intercept=FALSE)

    x.test.whole = cbind(1, X_test)
    # combine both S- and T- results
    best.beta = stats::coef(fit.tau,s="lambda.min")
    best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
    C.resST = x.test.whole %*% best.beta %*% W
    # based on S-
    best.beta = stats::coef(fit.tauS, s="lambda.min")
    best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
    C.resS = x.test.whole %*% best.beta %*% W
    # based on T-
    best.beta = stats::coef(fit.tauT, s="lambda.min")
    best.beta = matrix(best.beta[-1], nrow = pobs+1, byrow = T)
    C.resT = x.test.whole %*% best.beta %*% W
    # deC-learner results
    C.res = list(C.resST = C.resST, C.resS = C.resS, C.resT = C.resT)

    ######=================  General: R-learner  =====================######
    if ("R" %in% Learners) {
      # R will have multiple results with different reference group
      R.res = list()
      for (r in seq(length(K.grp))) {
        pi.hat = ps.train[, -r]
        W.ind = sapply(K.grp[-r], function(z){as.numeric(Trt==z)})
        adj.pi = W.ind - pi.hat
        x.tilde = matrix(as.vector(apply(adj.pi, 2, function(j) j * cbind(1,X_train))), nrow = nrow(X_train), byrow = F)
        penalty_factor = rep(c(0, array(1, pobs)), length(K.grp)-1)
        fit.R = glmnet::cv.glmnet(x.tilde, y.adj, penalty.factor = penalty_factor, family = "gaussian", parallel = TRUE, maxit = 100000, intercept = F)

        ## test data
        Est.R = NULL
        for ( k in 1:(length(K.grp)-1) ) {
          ind = rep(0, length(K.grp)-1); ind[k] = 1
          newX = t(ind) %x% cbind(1,X_test)
          est.R = predict(fit.R, newx = newX, type = "response", s = "lambda.min")
          Est.R = cbind(Est.R, est.R)
        }
        colnames(Est.R) <- paste(K.grp[-which(K.grp == r)], r, sep = "-")
        R.res = c(R.res, list(Est.R))
      }
      names(R.res) <- paste0("R", K.grp)
    }

    ######=================  General: Simplex R-learner (Rsim)  =====================######
    if ("Rsim" %in% Learners) {
      Rsim.res <- Rsim(X_train, y.adj, Trt, X_test, ps.train, ps.test, W)
    }


    ######=================  General: AD-learner  =====================######
    if ("AD" %in% Learners) {
      wts    = 1/sapply(seq(n.train), function(i) ps.train[i, Trt[i]])
      Y.mat  = sapply(seq(n.train), function(i) length(K.grp)*Y[i]*W[,Trt[i]])
      AD.fit = cv.glmnet(X_train, t(Y.mat), family="mgaussian", weights = wts, parallel = TRUE, intercept = TRUE)
      AD.res = predict(AD.fit, newx = X_test, s = "lambda.min")[,,1] %*% W
    }

    # returns
    colnames(S.res) <- colnames(T.res) <- colnames(Rsim.res) <- colnames(C.res$C.resST) <- colnames(C.res$C.resS) <- colnames(C.res$C.resT) <- colnames(AD.res) <- K.grp
    return(list(S.res = S.res, T.res = T.res, X.res = X.res,
                R.res = R.res, Rsim.res = Rsim.res,
                C.res = C.res, AD.res = AD.res))
  }


}


