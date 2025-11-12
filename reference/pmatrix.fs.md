# Transition probability matrix from a fully-parametric, time-inhomogeneous Markov multi-state model

The transition probability matrix for time-inhomogeneous Markov
multi-state models fitted to time-to-event data with
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
This has \\r,s\\ entry giving the probability that an individual is in
state \\s\\ at time \\t\\, given they are in state \\r\\ at time \\0\\.

## Usage

``` r
pmatrix.fs(
  x,
  trans = NULL,
  t = 1,
  newdata = NULL,
  condstates = NULL,
  ci = FALSE,
  tvar = "trans",
  sing.inf = 1e+10,
  B = 1000,
  cl = 0.95,
  tidy = FALSE,
  ...
)
```

## Arguments

- x:

  A model fitted with
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md)
  for the required form of the model and the data. Additionally, this
  must be a Markov / clock-forward model, but can be time-inhomogeneous.
  See the package vignette for further explanation.

  `x` can also be a list of models, with one component for each
  permitted transition, as illustrated in
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- trans:

  Matrix indicating allowed transitions. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- t:

  Time or vector of times to predict state occupancy probabilities for.

- newdata:

  A data frame specifying the values of covariates in the fitted model,
  other than the transition number. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- condstates:

  xInstead of the unconditional probability of being in state \\s\\ at
  time \\t\\ given state \\r\\ at time 0, return the probability
  conditional on being in a particular subset of states at time \\t\\.
  This subset is specified in the `condstates` argument, as a vector of
  character labels or integers.

  This is used, for example, in competing risks situations, e.g. if the
  competing states are death or recovery from a disease, and we want to
  compute the probability a patient has died, given they have died or
  recovered. If these are absorbing states, then as \\t\\ increases,
  this converges to the case fatality ratio. To compute this, set \\t\\
  to a very large number, `Inf` will not work.

- ci:

  Return a confidence interval calculated by simulating from the
  asymptotic normal distribution of the maximum likelihood estimates.
  Turned off by default, since this is computationally intensive. If
  turned on, users should increase `B` until the results reach the
  desired precision.

- tvar:

  Variable in the data representing the transition type. Not required if
  `x` is a list of models.

- sing.inf:

  If there is a singularity in the observed hazard, for example a
  Weibull distribution with `shape < 1` has infinite hazard at `t=0`,
  then as a workaround, the hazard is assumed to be a large finite
  number, `sing.inf`, at this time. The results should not be sensitive
  to the exact value assumed, but users should make sure by adjusting
  this parameter in these cases.

- B:

  Number of simulations from the normal asymptotic distribution used to
  calculate variances. Decrease for greater speed at the expense of
  accuracy.

- cl:

  Width of symmetric confidence intervals, relative to 1.

- tidy:

  If TRUE then return the results as a tidy data frame

- ...:

  Arguments passed to [`ode`](https://rdrr.io/pkg/deSolve/man/ode.html)
  in deSolve.

## Value

The transition probability matrix, if `t` is of length 1. If `t` is
longer, return a list of matrices, or a data frame if `tidy` is TRUE.

If `ci=TRUE`, each element has attributes `"lower"` and `"upper"` giving
matrices of the corresponding confidence limits. These are formatted for
printing but may be extracted using
[`attr()`](https://rdrr.io/r/base/attr.html).

## Details

This is computed by solving the Kolmogorov forward differential equation
numerically, using the methods in the deSolve package. The equation is

\$\$\frac{dP(t)}{dt} = P(t) Q(t)\$\$

where \\P(t)\\ is the transition probability matrix for time \\t\\, and
\\Q(t)\\ is the transition hazard or intensity as a function of \\t\\.
The initial condition is \\P(0) = I\\.

Note that the package msm has a similar method `pmatrix.msm`.
`pmatrix.fs` should give the same results as `pmatrix.msm` when both of
these conditions hold:

- the time-to-event distribution is exponential for all transitions,
  thus the `flexsurvreg` model was fitted with `dist="exp"` and the
  model is time-homogeneous.

- the msm model was fitted with `exacttimes=TRUE`, thus all the event
  times are known, and there are no time-dependent covariates.

msm only allows exponential or piecewise-exponential time-to-event
distributions, while flexsurvreg allows more flexible models. msm
however was designed in particular for panel data, where the process is
observed only at arbitrary times, thus the times of transition are
unknown, which makes flexible models difficult.

This function is only valid for Markov ("clock-forward") multi-state
models, though no warning or error is currently given if the model is
not Markov. See
[`pmatrix.simfs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.simfs.md)
for the equivalent for semi-Markov ("clock-reset") models.

## See also

[`pmatrix.simfs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.simfs.md),
[`totlos.fs`](http://chjackson.github.io/flexsurv-dev/reference/totlos.fs.md),
[`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.

## Examples

``` r
# BOS example in vignette, and in msfit.flexsurvreg
bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans,
                    data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
# more likely to be dead (state 3) as time moves on, or if start with
# BOS (state 2)
pmatrix.fs(bexp, t=c(5,10), trans=tmat)
#> $`5`
#>           [,1]      [,2]      [,3]
#> [1,] 0.2962297 0.2672185 0.4365518
#> [2,] 0.0000000 0.2672312 0.7327688
#> [3,] 0.0000000 0.0000000 1.0000000
#> 
#> $`10`
#>            [,1]       [,2]      [,3]
#> [1,] 0.08775208 0.15056691 0.7616810
#> [2,] 0.00000000 0.07141257 0.9285874
#> [3,] 0.00000000 0.00000000 1.0000000
#> 
#> attr(,"nst")
#> [1] 3
```
