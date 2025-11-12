# Transition probability matrix from a fully-parametric, semi-Markov multi-state model

The transition probability matrix for semi-Markov multi-state models
fitted to time-to-event data with
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
This has \\r,s\\ entry giving the probability that an individual is in
state \\s\\ at time \\t\\, given they are in state \\r\\ at time \\0\\.

## Usage

``` r
pmatrix.simfs(
  x,
  trans,
  t = 1,
  newdata = NULL,
  ci = FALSE,
  tvar = "trans",
  tcovs = NULL,
  M = 1e+05,
  B = 1000,
  cl = 0.95,
  cores = NULL,
  tidy = FALSE
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

  Time to predict state occupancy probabilities for. This can be a
  single number or a vector of different numbers.

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

- tidy:

  If `TRUE` then the results are returned as a tidy data frame with
  columns for the estimate and confidence limits, and rows per state
  transition and time interval.

## Value

The transition probability matrix. If `ci=TRUE`, there are attributes
`"lower"` and `"upper"` giving matrices of the corresponding confidence
limits. These are formatted for printing but may be extracted using
[`attr()`](https://rdrr.io/r/base/attr.html).

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
very inefficient, so that looping over the `M` individuals will be very
slow.

[`pmatrix.fs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.fs.md)
is a more efficient method based on solving the Kolmogorov forward
equation numerically, which requires the multi-state model to be Markov.
No error or warning is given if running `pmatrix.simfs` with a Markov
model, but this is still invalid.

## See also

[`pmatrix.fs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.fs.md),[`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md),[`totlos.simfs`](http://chjackson.github.io/flexsurv-dev/reference/totlos.simfs.md),
[`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.

## Examples

``` r
# BOS example in vignette, and in msfit.flexsurvreg

bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))

# more likely to be dead (state 3) as time moves on, or if start with
# BOS (state 2)

pmatrix.simfs(bexp, t=5, trans=tmat)
#>         1       2       3
#> 1 0.29456 0.26691 0.43853
#> 2 0.00000 0.26781 0.73219
#> 3 0.00000 0.00000 1.00000
pmatrix.simfs(bexp, t=10, trans=tmat)
#>         1       2       3
#> 1 0.08806 0.15161 0.76033
#> 2 0.00000 0.07006 0.92994
#> 3 0.00000 0.00000 1.00000

# these results should converge to those in help(pmatrix.fs), as M
# increases here and ODE solving precision increases there, since with
# an exponential distribution, the semi-Markov model is the same as the
# Markov model.
```
