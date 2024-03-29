% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetaLearners.R
\name{MetaLearners}
\alias{MetaLearners}
\title{Meta-learners for Treatment Recommendation and Individual Treatment Effect Estimation for multiple-treatment scenario}
\usage{
MetaLearners(
  X,
  Y,
  Trt,
  X.test = NULL,
  Learners = c("S", "T", "X", "R", "Rsim", "AD", "deC"),
  algorithm = "BART",
  controls = list(SL.library = c("SL.bartMachine", "SL.gam", "SL.ranger"), SL.library.PS
    = c("SL.bartMachine", "SL.gam", "SL.ranger"), basis.func = "Polynomial", degree = 2,
    n.knots = 3, X_train = NULL, X_test = NULL),
  ...
)
}
\arguments{
\item{X}{Data matrix}

\item{Y}{Outcome. Support either continuous or binary. If \code{Y} only has two levels, will be deemed
as binary outcome}

\item{Trt}{Treatment assignment indicator}

\item{X.test}{Test data. If \code{NULL}, training dataset, i.e., \code{X} will be used. In practice, \code{X.test} is always \code{NULL}.}

\item{Learners}{One or more of \code{"S", "T", "X", "R", "Rsim", "AD", "deC"}. Default includes all of them.
Notably, if "deC" is called, results from S- and T-learner will be reported as by product}

\item{algorithm}{Machine learning algorithm for the analysis. Currently, it supports one of \code{c("BART", "GAM", "RF", "SL")}.
representing Bayesian Additive Regression Trees (BART), Generalized Boosted Method (GBM), Random Forest (RF),
and Super Learner (SL). The default is \code{BART}.}

\item{controls}{A list of control arguments for analysis. See details.}

\item{...}{Additional argument for corresponding base-learners, e.g., \code{bart} in \code{dbarts} for BART,
\code{ranger} in package \code{ranger} for random forest.}
}
\value{
\item{S.res}{Results from S-learner. If binary outcome, then it represents the probability of \code{T=1}.
Causal estimands like log(Relative Risk), or log(Odds Ratio) can be obtained
from these estimated probabilities. Same for the rest.}
\item{T.res}{Results from T-learner}
\item{X.res}{Results from X-learner}
\item{R.res}{Results from R-learner. Notably, it contains multiple list of results with different
choice of reference. The results can be different with the reference level selection, especially
when treatment allocation is unbalanced.}
\item{Rsim.res}{Results from simplex R-learner}
\item{C.res}{Results from deC-learner, which contains \code{C.resS}, \code{C.resT}, and \code{C.resST},
representing using S-, T-, and average of S- and T- in calculating the model average. When sample size
is small or treatment allocation is unbalanced, \code{C.resS} is more preferred. Notably, if binary outcome,
the reported values \code{Pr[Y=1]} are not the causal estimand deC-learner is designed
for; please transform them to log(OR).}
\item{AD.res}{Results from AD-learning}
}
\description{
ITE/HTE/CATE estimation can be achieved by S-, T-, X-, R-, deC-learning,
while optimal treatment can be achieved by S-, T-, Rsim-, deC- and AD-learning.
}
\details{
For S-, T-, X-, and R-learner, they are initially proposed for two-treatment scenarios,
but can be extended to multiple-treatment following the same gist. R-learner will yield multiple sets
of results with the choise of reference group, and the results can be different from each other,
especially when treatment allocation is unbalanced. This algorithm will return results from all possible
reference group, named as \code{R.k}. User can then decide which reference group to take (technically, the group with largest
sample is more suitable as reference one). \cr\cr
R-, simplex R (Rsim), deC-, and AD-learning have to choose the basis function since a global additive structure
is assumed. Different basis functions are available and user can adjust the corresponding \code{controls} arguments.
Interaction terms are not available to generate automatically but user can input their own design matrix.
The lasso-type penalization from \code{glmnet} will be adopted in analysis for R-, simplex R (Rsim), and deC-learner,
and group lasso for AD-learning.
\cr\cr
\code{controls} is a list containing:
\itemize{
\item{SL.library}{: The algorithms to be included in for Super Learner. For algorithm is \code{SL} only.
For details, please check package \link{SuperLearner}}
\item{SL.library.PS}{: The SL.library for propensity score estimation. It can be different from \code{SL.library}
given some algorithms may not capable for multinomial outcomes. For algorithm is \code{SL} only.}
\item{basis.func}{: To specify the basis functions. Currently support one of these
\code{c("Polynomial", "bs", "ns")}, standing for polynomial regression, B-splines, and natural cubic spline.}
\item{degree}{: The degree for basis functions (not for natural cubic spline)}
\item{n.knots}{: Number of knots when B-slines or natural spline is adopted.}
\item{X_train}{: User provided design matrix. If given, it overrides all basis function related inputs.}
\item{X_test}{: Similar to \code{X_train}. For really world application rather than simulation studies,
it is not needed.}
}
}
\examples{
\dontrun{
require(doMC)
registerDoMC(cores = 6)
res = MetaLearners(simu.dat$X, simu.dat$Y, simu.dat$Trt, algorithm = "GAM")
final.res = parseRes(res)
}
}
