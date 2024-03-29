# Meta-Learners
The package is developed for treatment recommendation &amp; pairwise treatment individual effect estimation (ITE/CATE/HTE). It includes some published methods such as S-learner, T-learner, X-learner, and R-learner and some newer/unpublished methods like reference-free simplex R-learner and de-Centralized-learner. Angle-based direct learning (AD-learner), an outcome-weighting-based treatment recommendation method is included as well. 

Notably, the S-, T-, and X-learner follows the nomenclature of [Kunzel et. al](https://www.pnas.org/doi/pdf/10.1073/pnas.1804597116). There are many causal inference methods follows the S- and T-learner structure, like causal boosting, causal forests, etc. X-learner is designed for two-treatment setting and is extended for multiple treatment case in the package.

R-learner is proposed by [Nie and Wager](https://academic.oup.com/biomet/article/108/2/299/5911092?login=true) mainly for two-treatment setting. Here we extend it for multiple treatment setting as well but we found that different reference group may cause different results. 

AD-learning is relatively new method by [Qi et. al.](https://www.tandfonline.com/doi/abs/10.1080/01621459.2018.1529597) which does not involve the estimation of treatment effects but ITR. Individualized Treatment Rules/Regimens/Recommendations (ITR), a mapping from covariate space to treatment space, is a treatment decision rule that determines the optimal treatment directly given the subject's covariates. 

Reference-free R-learner is proposed by [Zhou et. al.](https://journals.sagepub.com/doi/10.1177/09622802221144326) which follows the idea of R-learner but fixes the issue of recommendation inconsistency in multi-armed scenarios.

de-Centralized-learner is now under development. For more details, please wait for publication.

Some technical details can be found on author's [website](https://jzhou.org/posts/metalearner/).

# Installation
To install the package:
```
devtools::install_github("junyzhou10/MetaLearners")
```

# Updates (11/13/2022)
Now binary outcomes is allowed for S-, T-, and de-Centralized learner. Notice that, binary outcome is not supposed to be supported by X-, and R-learner methods. AD-learner is capable, and will be developed soon.
