# Simulate paths through a fully parametric semi-Markov multi-state model

Simulate changes of state and transition times from a semi-Markov
multi-state model fitted using
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

## Usage

``` r
sim.fmsm(
  x,
  trans = NULL,
  t,
  newdata = NULL,
  start = 1,
  M = 10,
  tvar = "trans",
  tcovs = NULL,
  tidy = FALSE
)
```

## Arguments

- x:

  A model fitted with
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md)
  for the required form of the model and the data.

  Alternatively `x` can be a list of fitted
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  model objects. The `i`th element of this list is the model
  corresponding to the `i`th transition in `trans`. This is a more
  efficient way to fit a multi-state model, but only valid if the
  parameters are different between different transitions.

- trans:

  Matrix indicating allowed transitions. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- t:

  Time, or vector of times for each of the `M` individuals, to simulate
  trajectories until.

- newdata:

  A data frame specifying the values of covariates in the fitted model,
  other than the transition number. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- start:

  Starting state, or vector of starting states for each of the `M`
  individuals.

- M:

  Number of individual trajectories to simulate.

- tvar:

  Variable in the data representing the transition type. Not required if
  `x` is a list of models.

- tcovs:

  Names of "predictable" time-dependent covariates in `newdata`, i.e.
  those whose values change at the same rate as time. Age is a typical
  example. During simulation, their values will be updated after each
  transition time, by adding the current time to the value supplied in
  `newdata`. This assumes the covariate is measured in the same unit as
  time. `tcovs` is supplied as a character vector.

- tidy:

  If `TRUE` then the simulated data are returned as a tidy data frame
  with one row per simulated transition. See
  [`simfs_bytrans`](http://chjackson.github.io/flexsurv-dev/reference/simfs_bytrans.md)
  for a function to rearrange the data into this format if it was
  simulated in non-tidy format. Currently this includes only event
  times, and excludes any times of censoring that are reported when
  `tidy=FALSE`.

## Value

If `tidy=TRUE`, a data frame with one row for each simulated transition,
giving the individual ID `id`, start state `start`, end state `end`,
transition label `trans`, time of the transition since the start of the
process (`time`), and time since the previous transition (`delay`).

If `tidy=FALSE`, a list of two matrices named `st` and `t`. The rows of
each matrix represent simulated individuals. The columns of `t` contain
the times when the individual changes state, to the corresponding states
in `st`.

The first columns will always contain the starting states and the
starting times. The last column of `t` represents either the time when
the individual moves to an absorbing state, or right-censoring in a
transient state at the time given in the `t` argument to `sim.fmsm`.

## Details

`sim.fmsm` relies on the presence of a function to sample random numbers
from the parametric survival distribution used in the fitted model `x`,
for example [`rweibull`](https://rdrr.io/r/stats/Weibull.html) for
Weibull models. If `x` was fitted using a custom distribution, called
`dist` say, then there must be a function called (something like)
`rdist` either in the working environment, or supplied through the
`dfns` argument to
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
This must be in the same format as standard R functions such as
[`rweibull`](https://rdrr.io/r/stats/Weibull.html), with first argument
`n`, and remaining arguments giving the parameters of the distribution.
It must be vectorised with respect to the parameter arguments.

This function is only valid for semi-Markov ("clock-reset") models,
though no warning or error is currently given if the model is not of
this type. An equivalent for time-inhomogeneous Markov ("clock-forward")
models has currently not been implemented.

## See also

[`pmatrix.simfs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.simfs.md),[`totlos.simfs`](http://chjackson.github.io/flexsurv-dev/reference/totlos.simfs.md)

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.

## Examples

``` r
bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
sim.fmsm(bexp, M=10, t=5, trans=tmat)
#> $st
#>       [,1] [,2] [,3]
#>  [1,]    1    2    3
#>  [2,]    1    2    3
#>  [3,]    1    2    2
#>  [4,]    1    2    3
#>  [5,]    1    2    2
#>  [6,]    1    2    3
#>  [7,]    1    2    2
#>  [8,]    1    2    2
#>  [9,]    1    2    2
#> [10,]    1    2    2
#> 
#> $t
#>       [,1]       [,2]      [,3]
#>  [1,]    0 0.05740866 0.2204663
#>  [2,]    0 2.77917357 3.4298756
#>  [3,]    0 3.12006461 5.0000000
#>  [4,]    0 2.26788620 2.5892033
#>  [5,]    0 0.59817485 5.0000000
#>  [6,]    0 0.10337274 1.3235881
#>  [7,]    0 3.04681142 5.0000000
#>  [8,]    0 3.33806361 5.0000000
#>  [9,]    0 2.74151326 5.0000000
#> [10,]    0 1.19843599 5.0000000
#> 
#> attr(,"trans")
#>      [,1] [,2] [,3]
#> [1,]   NA    1    2
#> [2,]   NA   NA    3
#> [3,]   NA   NA   NA
```
