# Expected total length of stay in specific states, from a fully-parametric, semi-Markov multi-state model

The expected total time spent in each state for semi-Markov multi-state
models fitted to time-to-event data with
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
This is defined by the integral of the transition probability matrix,
though this is not analytically possible and is computed by simulation.

## Usage

``` r
totlos.simfs(
  x,
  trans,
  t = 1,
  start = 1,
  newdata = NULL,
  ci = FALSE,
  tvar = "trans",
  tcovs = NULL,
  group = NULL,
  M = 1e+05,
  B = 1000,
  cl = 0.95,
  cores = NULL
)
```

## Arguments

- x:

  A model fitted with
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md)
  for the required form of the model and the data. Additionally this
  should be semi-Markov, so that the time variable represents the time
  since the last transition. In other words the response should be of
  the form `Surv(time,status)`. See the package vignette for further
  explanation.

  `x` can also be a list of
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  models, with one component for each permitted transition, as
  illustrated in
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).
  This can be constructed by
  [`fmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmsm.md).

- trans:

  Matrix indicating allowed transitions. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).
  This is not required if `x` is a list constructed by
  [`fmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmsm.md).

- t:

  Maximum time to predict to.

- start:

  Starting state.

- newdata:

  A data frame specifying the values of covariates in the fitted model,
  other than the transition number. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- ci:

  Return a confidence interval calculated by simulating from the
  asymptotic normal distribution of the maximum likelihood estimates.
  This is turned off by default, since two levels of simulation are
  required. If turned on, users should adjust `B` and/or `M` until the
  results reach the desired precision. The simulation over `M` is
  generally vectorised, therefore increasing `B` is usually more
  expensive than increasing `M`.

- tvar:

  Variable in the data representing the transition type. Not required if
  `x` is a list of models.

- tcovs:

  Predictable time-dependent covariates such as age, see
  [`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md).

- group:

  Optional grouping for the states. For example, if there are four
  states, and `group=c(1,1,2,2)`, then `totlos.simfs` returns the
  expected total time in states 1 and 2 combined, and states 3 and 4
  combined.

- M:

  Number of individuals to simulate in order to approximate the
  transition probabilities. Users should adjust this to obtain the
  required precision.

- B:

  Number of simulations from the normal asymptotic distribution used to
  calculate confidence limits. Decrease for greater speed at the expense
  of accuracy.

- cl:

  Width of symmetric confidence intervals, relative to 1.

- cores:

  Number of processor cores used when calculating confidence limits by
  repeated simulation. The default uses single-core processing.

## Value

The expected total time spent in each state (or group of states given by
`group`) up to time `t`, and corresponding confidence intervals if
requested.

## Details

This is computed by simulating a large number of individuals `M` using
the maximum likelihood estimates of the fitted model and the function
[`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md).
Therefore this requires a random sampling function for the parametric
survival model to be available: see the "Details" section of
[`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md).
This will be available for all built-in distributions, though users may
need to write this for custom models.

Note the random sampling method for `flexsurvspline` models is currently
very inefficient, so that looping over `M` will be very slow.

The equivalent function for time-inhomogeneous Markov models is
[`totlos.fs`](http://chjackson.github.io/flexsurv-dev/reference/totlos.fs.md).
Note neither of these functions give errors or warnings if used with the
wrong type of model, but the results will be invalid.

## See also

[`pmatrix.simfs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.simfs.md),[`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md),[`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.

## Examples

``` r
# BOS example in vignette, and in msfit.flexsurvreg
bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))

# predict 4 years spent without BOS, 3 years with BOS, before death
# As t increases, this should converge
totlos.simfs(bexp, t=10, trans=tmat)
#>        1        2        3 
#> 3.744624 2.124969 4.130406 
totlos.simfs(bexp, t=1000, trans=tmat)
#>          1          2          3 
#>   4.127459   2.957657 992.914884 
```
