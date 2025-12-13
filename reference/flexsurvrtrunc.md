# Flexible parametric models for right-truncated, uncensored data defined by times of initial and final events.

This function estimates the distribution of the time between an initial
and final event, in situations where individuals are only observed if
they have experienced both events before a certain time, thus they are
right-truncated at this time. The time of the initial event provides
information about the time from initial to final event, given the
truncated observation scheme, and initial events are assumed to occur
with an exponential growth rate.

## Usage

``` r
flexsurvrtrunc(
  t,
  tinit,
  rtrunc,
  tmax,
  data = NULL,
  method = "joint",
  dist,
  theta = NULL,
  fixed.theta = TRUE,
  inits = NULL,
  fixedpars = NULL,
  dfns = NULL,
  integ.opts = NULL,
  cl = 0.95,
  optim_control = list()
)
```

## Arguments

- t:

  Vector of time differences between an initial and final event for a
  set of individuals.

- tinit:

  Absolute time of the initial event for each individual.

- rtrunc:

  Individual-specific right truncation points on the same scale as `t`,
  so that each individual's survival time `t` would not have been
  observed if it was greater than the corresponding element of `rtrunc`.
  Only used in `method="joint"`. In `method="final"`, the
  right-truncation is implicit.

- tmax:

  Maximum possible time between initial and final events that could have
  been observed. This is only used in `method="joint"`, and is ignored
  in `method="final"`.

- data:

  Data frame containing `t`, `rtrunc` and `tinit`.

- method:

  If `"joint"` then the "joint-conditional" method is used. If `"final"`
  then the "conditional-on-final" method is used. The
  "conditional-on-initial" method can be implemented by using
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  with a `rtrunc` argument. These methods are all described in Seaman et
  al. (2020).

- dist:

  Typically, one of the strings in the first column of the following
  table, identifying a built-in distribution. This table also identifies
  the location parameters, and whether covariates on these parameters
  represent a proportional hazards (PH) or accelerated failure time
  (AFT) model. In an accelerated failure time model, the covariate
  speeds up or slows down the passage of time. The relation of the
  reported coefficient to the "time acceleration factor" is described
  for each distribution in the "Distributions reference" vignette. Note
  the interpretation of the sign of the coefficient may be different for
  different distributions.

  |                   |                              |         |     |
  |-------------------|------------------------------|---------|-----|
  | `"gengamma"`      | Generalized gamma (stable)   | mu      | AFT |
  | `"gengamma.orig"` | Generalized gamma (original) | scale   | AFT |
  | `"genf"`          | Generalized F (stable)       | mu      | AFT |
  | `"genf.orig"`     | Generalized F (original)     | mu      | AFT |
  | `"weibull"`       | Weibull                      | scale   | AFT |
  | `"weibullPH"`     | Weibull                      | scale   | PH  |
  | `"gamma"`         | Gamma                        | rate    | AFT |
  | `"exp"`           | Exponential                  | rate    | PH  |
  | `"llogis"`        | Log-logistic                 | scale   | AFT |
  | `"lnorm"`         | Log-normal                   | meanlog | AFT |
  | `"gompertz"`      | Gompertz                     | rate    | PH  |

  `"exponential"` and `"lognormal"` can be used as aliases for `"exp"`
  and `"lnorm"`, for compatibility with
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

  Alternatively, `dist` can be a list specifying a custom distribution.
  See section “Custom distributions” below for how to construct this
  list.

  Very flexible spline-based distributions can also be fitted with
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

  The parameterisations of the built-in distributions used here are the
  same as in their built-in distribution functions:
  [`dgengamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md),
  [`dgengamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md),
  [`dgenf`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md),
  [`dgenf.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenF.orig.md),
  [`dweibull`](https://rdrr.io/r/stats/Weibull.html),
  [`dgamma`](https://rdrr.io/r/stats/GammaDist.html),
  [`dexp`](https://rdrr.io/r/stats/Exponential.html),
  [`dlnorm`](https://rdrr.io/pkg/eha/man/Lognormal.html),
  [`dgompertz`](http://chjackson.github.io/flexsurv-dev/reference/Gompertz.md),
  respectively. The functions in base R are used where available,
  otherwise, they are provided in this package.

  A package vignette "Distributions reference" lists the survivor
  functions and covariate effect parameterisations used by each built-in
  distribution.

  For the Weibull, exponential and log-normal distributions,
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  simply works by calling
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) to obtain
  the maximum likelihood estimates, then calling
  [`optim`](https://rdrr.io/r/stats/optim.html) to double-check
  convergence and obtain the covariance matrix for
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)'s
  preferred parameterisation.

  The Weibull parameterisation is different from that in
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html), instead it
  is consistent with [`dweibull`](https://rdrr.io/r/stats/Weibull.html).
  The `"scale"` reported by
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) is
  equivalent to `1/shape` as defined by
  [`dweibull`](https://rdrr.io/r/stats/Weibull.html) and hence
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  The first coefficient `(Intercept)` reported by
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) is
  equivalent to `log(scale)` in
  [`dweibull`](https://rdrr.io/r/stats/Weibull.html) and
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

  Similarly in the exponential distribution, the rate, rather than the
  mean, is modelled on covariates.

  The object `flexsurv.dists` lists the names of the built-in
  distributions, their parameters, location parameter, functions used to
  transform the parameter ranges to and from the real line, and the
  functions used to generate initial values of each parameter for
  estimation.

- theta:

  Initial value (or fixed value) for the exponential growth rate
  `theta`. Defaults to 1.

- fixed.theta:

  Should `theta` be fixed at its initial value or estimated. This only
  applies to `method="joint"`. For `method="final"`, `theta` must be
  fixed.

- inits:

  Initial values for the parameters of the parametric survival
  distributon. If not supplied, a heuristic is used. as is done in
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- fixedpars:

  Integer indices of the parameters of the survival distribution that
  should be fixed to their values supplied in `inits`. Same length as
  `inits`.

- dfns:

  An alternative way to define a custom survival distribution (see
  section “Custom distributions” below). A list whose components may
  include `"d"`, `"p"`, `"h"`, or `"H"` containing the probability
  density, cumulative distribution, hazard, or cumulative hazard
  functions of the distribution. For example,

  `list(d=dllogis, p=pllogis)`.

  If `dfns` is used, a custom `dlist` must still be provided, but
  `dllogis` and `pllogis` need not be visible from the global
  environment. This is useful if `flexsurvreg` is called within other
  functions or environments where the distribution functions are also
  defined dynamically.

- integ.opts:

  List of named arguments to pass to
  [`integrate`](https://rdrr.io/r/stats/integrate.html), if a custom
  density or hazard is provided without its cumulative version. For
  example,

  `integ.opts = list(rel.tol=1e-12)`

- cl:

  Width of symmetric confidence intervals for maximum likelihood
  estimates, by default 0.95.

- optim_control:

  List to supply as the `control` argument to
  [`optim`](https://rdrr.io/r/stats/optim.html) to control the
  likelihood maximisation.

## Details

Covariates are not currently supported.

Note that
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
with an `rtrunc` argument, can fit models for a similar kind of data,
but those models ignore the information provided by the time of the
initial event.

A nonparametric estimator of survival under right-truncation is also
provided in
[`survrtrunc`](http://chjackson.github.io/flexsurv-dev/reference/survrtrunc.md).
See Seaman et al. (2020) for a full comparison of the alternative
models.

## References

Seaman, S., Presanis, A. and Jackson, C. (2020) Estimating a
Time-to-Event Distribution from Right-Truncated Data in an Epidemic: a
Review of Methods

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`survrtrunc`](http://chjackson.github.io/flexsurv-dev/reference/survrtrunc.md).

## Examples

``` r
set.seed(1) 
## simulate time to initial event
X <- rexp(1000, 0.2)
## simulate time between initial and final event
T <- rgamma(1000, 2, 10) 

## right-truncate to keep only those with final event time
## before tmax
tmax <- 40
obs <- X + T < tmax 
rtrunc <- tmax - X
dat <- data.frame(X, T, rtrunc)[obs,]

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gamma", theta=0.2)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "gamma", theta = 0.2)
#> 
#> Estimates: 
#>         pars   est   lcl   ucl  estlog   selog
#> shape  shape  1.93  1.78   2.1   0.659  0.0414
#> rate    rate  9.22  8.39  10.1   2.222  0.0483
#> theta  theta  1.22    NA    NA   0.200      NA
#> 
#> Log-likelihood = -7846.25, df = 2
#> AIC = 15696.5
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gamma", theta=0.2, fixed.theta=FALSE)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "gamma", theta = 0.2, fixed.theta = FALSE)
#> 
#> Estimates: 
#>         pars    est    lcl     ucl  estlog    selog
#> shape  shape  1.932  1.781   2.095   0.658  0.04144
#> rate    rate  9.419  8.586  10.334   2.243  0.04728
#> theta  theta  0.824  0.814   0.834  -0.193  0.00621
#> 
#> Log-likelihood = -1949.284, df = 3
#> AIC = 3904.568
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gamma", theta=0.2, inits=c(1, 8))
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "gamma", theta = 0.2, inits = c(1, 8))
#> 
#> Estimates: 
#>         pars   est   lcl   ucl  estlog   selog
#> shape  shape  1.93  1.78   2.1   0.659  0.0414
#> rate    rate  9.22  8.39  10.1   2.222  0.0483
#> theta  theta  1.22    NA    NA   0.200      NA
#> 
#> Log-likelihood = -7846.25, df = 2
#> AIC = 15696.5
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gamma", theta=0.2, method="final")
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, method = "final", dist = "gamma", theta = 0.2)
#> 
#> Estimates: 
#>         pars   est   lcl   ucl  estlog   selog
#> shape  shape  1.94  1.79  2.11   0.663  0.0415
#> rate    rate  9.05  8.22  9.97   2.203  0.0493
#> theta  theta  1.22    NA    NA   0.200      NA
#> 
#> Log-likelihood = 715.5208, df = 2
#> AIC = -1427.042
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gamma", fixed.theta=TRUE)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "gamma", fixed.theta = TRUE)
#> 
#> Estimates: 
#>         pars   est   lcl   ucl  estlog   selog
#> shape  shape  1.93  1.78  2.09   0.658  0.0415
#> rate    rate  8.41  7.58  9.33   2.130  0.0529
#> theta  theta  2.72    NA    NA   1.000      NA
#> 
#> Log-likelihood = -33947.91, df = 2
#> AIC = 67899.82
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="weibull", fixed.theta=TRUE)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "weibull", fixed.theta = TRUE)
#> 
#> Estimates: 
#>         pars    est    lcl    ucl  estlog   selog
#> shape  shape  1.501  1.430  1.577   0.406  0.0249
#> scale  scale  0.252  0.241  0.264  -1.378  0.0240
#> theta  theta  2.718     NA     NA   1.000      NA
#> 
#> Log-likelihood = -33950.21, df = 2
#> AIC = 67904.42
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="lnorm", fixed.theta=TRUE)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "lnorm", fixed.theta = TRUE)
#> 
#> Estimates: 
#>             pars     est     lcl    ucl  estlog   selog
#> meanlog  meanlog  -1.694  -1.758  -1.63  -1.694  0.0325
#> sdlog      sdlog   0.893   0.848   0.94  -0.113  0.0264
#> theta      theta   2.718      NA     NA   1.000      NA
#> 
#> Log-likelihood = -33991.77, df = 2
#> AIC = 67987.55
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gengamma", fixed.theta=TRUE)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "gengamma", fixed.theta = TRUE)
#> 
#> Estimates: 
#>         pars     est     lcl     ucl  estlog   selog
#> mu        mu  -1.446  -1.520  -1.371  -1.446  0.0379
#> sigma  sigma   0.702   0.657   0.750  -0.354  0.0335
#> Q          Q   0.798   0.631   0.964   0.798  0.0851
#> theta  theta   2.718      NA      NA   1.000      NA
#> 
#> Log-likelihood = -33947.48, df = 3
#> AIC = 67900.97
#> 

flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                dist="gompertz", fixed.theta=TRUE)
#> Call:
#> flexsurvrtrunc(t = T, tinit = X, rtrunc = rtrunc, tmax = 40, 
#>     data = dat, dist = "gompertz", fixed.theta = TRUE)
#> 
#> Estimates: 
#>         pars   est   lcl   ucl  estlog   selog
#> shape  shape  2.58  2.18  2.97   2.576  0.2009
#> rate    rate  2.63  2.37  2.92   0.968  0.0532
#> theta  theta  2.72    NA    NA   1.000      NA
#> 
#> Log-likelihood = -33990.12, df = 2
#> AIC = 67984.25
#> 
```
