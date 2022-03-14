#' @title Parse results from MetaLearners()
#' @description Find the optimal treatment & pairwise treatment effect given results from MetaLearner()
#' @param Meta.Res Return object from [MetaLearners()]
#' @details Since results from [MetaLearners()] only provides the marginal outcomes for many learners,
#' so it is necessary to translate the results into pairwise CATE/HTE/ITE and optimal treatment recommendations.
#' Note that for X-learner, we only return the estimated pairwise treatment effect but no optimal results.
#' People can choose their own method to figure out the optimal one based on pairwise results.
#' For simplex R-learner (Rsim) and AD-learning, only optimal treatment is provided as they are not designed to
#' estimate pairwise treatment effect.
#'
#' @return A list of results where each element is a list of \code{ITE}, the estimated pairwise
#' treatment effect, and \code{Opt.trt}, the suggested optimal treatment.
#' @md
#' @export

parseRes <- function(Meta.Res) {
  # utility function
  util.func <- function(mat) {
    name.seq = colnames(mat)
    p = ncol(mat)
    res = NULL
    res.name = NULL
    for (i in seq(p)) {
      for (j in seq(p)[-(1:i)]) {
        res = cbind(res, mat[,i]-mat[,j])
        res.name = c(res.name, paste(name.seq[i],name.seq[j],sep="-"))
      }
    }
    colnames(res) <- res.name
    return(res)
  }
  # parse R-learner
  parseR <- function(mat) {
    name.seq = colnames(mat)
    name.seq = unlist(strsplit(name.seq, "-"))
    ref  = name.seq[2]
    trts = name.seq[seq(from = 1, to = length(name.seq), by = 2)]
    opt.trt = apply(mat, 1, function(x) {
      if (sum(x<0)==length(trts)) {
        return(ref)
      } else {
        return(trts[which.max(x)])
      }
    })
    return(opt.trt)
  }

  S.res <- T.res <- X.res <- R1.res <- R2.res <- R3.res <- R4.res <- Rsim.res <- C.resT <- C.resS <- C.resST <- AD.res <- list()
  name.seq = colnames(Meta.Res$S.res)

  # X-learner(no need to do adjust)
  X.res$ITE <- Meta.Res$X.res
  # S-learner
  S.res$Opt.trt = name.seq[apply(Meta.Res$S.res, 1, which.max)]
  S.res$ITE = util.func(Meta.Res$S.res)
  # T-learner
  T.res$Opt.trt = name.seq[apply(Meta.Res$T.res, 1, which.max)]
  T.res$ITE = util.func(Meta.Res$T.res)
  # deC-learners
  C.resT$Opt.trt = name.seq[apply(Meta.Res$C.res$C.resT, 1, which.max)]
  C.resT$ITE = util.func(Meta.Res$C.res$C.resT)
  C.resS$Opt.trt = name.seq[apply(Meta.Res$C.res$C.resS, 1, which.max)]
  C.resS$ITE = util.func(Meta.Res$C.res$C.resS)
  C.resST$Opt.trt = name.seq[apply(Meta.Res$C.res$C.resST, 1, which.max)]
  C.resST$ITE = util.func(Meta.Res$C.res$C.resST)

  # R-learners
  if (is.null(Meta.Res$R.res)) {
    R1.res <- R2.res <- R3.res <- R4.res <- NULL
  } else {
    R1.res$Opt.trt = parseR(Meta.Res$R.res$R1)
    R1.res$ITE = Meta.Res$R.res$R1
    R2.res$Opt.trt = parseR(Meta.Res$R.res$R2)
    R2.res$ITE = Meta.Res$R.res$R2
    R3.res$Opt.trt = parseR(Meta.Res$R.res$R3)
    R3.res$ITE = Meta.Res$R.res$R3
    R4.res$Opt.trt = parseR(Meta.Res$R.res$R4)
    R4.res$ITE = Meta.Res$R.res$R4
  }

  # Rsim
  if (is.null(Meta.Res$Rsim.res)) {
    Rsim.res = NULL
  } else {
    Rsim.res$Opt.trt = name.seq[apply(Meta.Res$Rsim.res, 1, which.max)]
  }

  # AD-learning
  if (is.null(Meta.Res$AD.res)) {
    AD.res = NULL
  } else {
    AD.res$Opt.trt = name.seq[apply(Meta.Res$AD.res, 1, which.max)]
  }

  return(list(S.res = S.res, T.res = T.res, X.res = X.res,
              R1.res = R1.res, R2.res = R2.res, R3.res = R3.res, R4.res = R4.res,
              Rsim.res = Rsim.res,
              C.resT = C.resT, C.resS = C.resS, C.resST = C.resST,
              AD.res = AD.res))
}


