flexsurv
========

The [flexsurv](http://cran.r-project.org/package=flexsurv) R package for parametric survival and multi-state modelling.


## Key features

* Parametric models for time-to-event (survival) data.  Data may be right-censored, and/or left-censored, and/or left-truncated.  Several built-in parametric distributions are available, including a very flexible model based on splines (Royston and Parmar). 

* Any user-defined parametric model can be employed by supplying R functions defining the distribution.

* Covariates can be included using a (log-)linear model on any parameter of any distribution.  This typically defines an accelerated failure time or proportional hazards model.

* Multi-state models for continuously-observed data can be defined by piecing together transition-specific parametric models of any kind.   (For intermittently-observed data, see instead the [msm](http://CRAN.R-project.org/package=msm) package.)


## Learn more 

[User guide (PDF)](https://chjackson.github.io/flexsurv/articles/flexsurv.pdf).

[Guide to multi-state modelling in flexsurv](https://chjackson.github.io/flexsurv/articles/multistate.pdf).

[Full reference manual](https://chjackson.github.io/flexsurv/reference/index.html) for all the package's functions.

[Paper in Journal of Statistical Software (Jackson 2016)](https://www.jstatsoft.org/article/view/v070i08).  Mostly the same as the user guide, but not kept up to date.


## Installation (stable CRAN version)
```r
install.packages("flexsurv")
```

## Installation (development version)

```r
install.packages("devtools") # if devtools not already installed
devtools::install_github('chjackson/flexsurv-dev')
```

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/chjackson/flexsurv-dev/workflows/R-CMD-check/badge.svg)](https://github.com/chjackson/flexsurv-dev/actions)
[![Codecov test coverage](https://codecov.io/gh/chjackson/flexsurv/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chjackson/flexsurv?branch=master)
<!-- badges: end -->
