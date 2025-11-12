# Total length of stay in particular states for a fully-parametric, time-inhomogeneous Markov multi-state model

The matrix whose \\r,s\\ entry is the expected amount of time spent in
state \\s\\ for a time-inhomogeneous, continuous-time Markov multi-state
process that starts in state \\r\\, up to a maximum time \\t\\. This is
defined as the integral of the corresponding transition probability up
to that time.

## Usage

``` r
totlos.fs(
  x,
  trans = NULL,
  t = 1,
  newdata = NULL,
  ci = FALSE,
  tvar = "trans",
  sing.inf = 1e+10,
  B = 1000,
  cl = 0.95,
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

  Time or vector of times to predict up to. Must be finite.

- newdata:

  A data frame specifying the values of covariates in the fitted model,
  other than the transition number. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

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

- ...:

  Arguments passed to [`ode`](https://rdrr.io/pkg/deSolve/man/ode.html)
  in deSolve.

## Value

The matrix of lengths of stay \\T(t)\\, if `t` is of length 1, or a list
of matrices if `t` is longer.

If `ci=TRUE`, each element has attributes `"lower"` and `"upper"` giving
matrices of the corresponding confidence limits. These are formatted for
printing but may be extracted using
[`attr()`](https://rdrr.io/r/base/attr.html).

The result also has an attribute `P` giving the transition probability
matrices, since these are unavoidably computed as a side effect. These
are suppressed for printing, but can be extracted with `attr(...,"P")`.

## Details

This is computed by solving a second order extension of the Kolmogorov
forward differential equation numerically, using the methods in the
deSolve package. The equation is expressed as a linear system

\$\$\frac{dT(t)}{dt} = P(t)\$\$ \$\$\frac{dP(t)}{dt} = P(t) Q(t)\$\$

and solved for \\T(t)\\ and \\P(t)\\ simultaneously, where \\T(t)\\ is
the matrix of total lengths of stay, \\P(t)\\ is the transition
probability matrix for time \\t\\, and \\Q(t)\\ is the transition hazard
or intensity as a function of \\t\\. The initial conditions are \\T(0) =
0\\ and \\P(0) = I\\.

Note that the package msm has a similar method `totlos.msm`. `totlos.fs`
should give the same results as `totlos.msm` when both of these
conditions hold:

- the time-to-event distribution is exponential for all transitions,
  thus the `flexsurvreg` model was fitted with `dist="exp"`, and is
  time-homogeneous.

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
[`totlos.simfs`](http://chjackson.github.io/flexsurv-dev/reference/totlos.simfs.md)
for the equivalent for semi-Markov ("clock-reset") models.

## See also

[`totlos.simfs`](http://chjackson.github.io/flexsurv-dev/reference/totlos.simfs.md),
[`pmatrix.fs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.fs.md),
[`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.

## Examples

``` r
# BOS example in vignette, and in msfit.flexsurvreg
bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans,
                    data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))

# predict 4 years spent without BOS, 3 years with BOS, before death
# As t increases, this should converge

totlos.fs(bexp, t=10, trans=tmat)
#>          [,1]     [,2]      [,3]
#> [1,] 3.749105 2.126570  4.124326
#> [2,] 0.000000 3.518329  6.481671
#> [3,] 0.000000 0.000000 10.000000
totlos.fs(bexp, t=1000, trans=tmat)
#>          [,1]     [,2]      [,3]
#> [1,] 4.109742 2.956493  992.9338
#> [2,] 0.000000 3.788904  996.2111
#> [3,] 0.000000 0.000000 1000.0000
totlos.fs(bexp, t=c(5,10), trans=tmat)
#> $`5`
#>          [,1]     [,2]     [,3]
#> [1,] 2.892316 1.068225 1.039459
#> [2,] 0.000000 2.776392 2.223608
#> [3,] 0.000000 0.000000 5.000000
#> 
#> $`10`
#>          [,1]     [,2]      [,3]
#> [1,] 3.749105 2.126570  4.124326
#> [2,] 0.000000 3.518329  6.481671
#> [3,] 0.000000 0.000000 10.000000
#> 

# Answers should match results in help(totlos.simfs) up to Monte Carlo
# error there / ODE solving precision here, since with an exponential
# distribution, the "semi-Markov" model there is the same as the Markov
# model here
```
