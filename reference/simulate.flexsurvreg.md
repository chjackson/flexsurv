# Simulate censored time-to-event data from a fitted flexsurvreg model

Simulate censored time-to-event data from a fitted flexsurvreg model

## Usage

``` r
# S3 method for class 'flexsurvreg'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  newdata = NULL,
  start = NULL,
  censtime = NULL,
  tidy = FALSE,
  ...
)
```

## Arguments

- object:

  Object returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- nsim:

  Number of simulations per row in `newdata`.

- seed:

  Random number seed. This is returned with the result of this function,
  as described in [`simulate`](https://rdrr.io/r/stats/simulate.html)
  for the `lm` method.

- newdata:

  Data frame defining alternative sets of covariate values to simulate
  with. If omitted, this defaults to the data originally used to fit the
  model.

- start:

  Optional left-truncation time or times. The returned survival, hazard
  or cumulative hazard will be conditioned on survival up to this time.
  Predicted times returned with `"rmst"`, `"mean"`, `"median"` or
  `"quantile"` will be times since time zero, not times since the
  `start` time.

  A vector of the same length as `t` can be supplied to allow different
  truncation times for each prediction time, though this doesn't make
  sense in the usual case where this function is used to calculate a
  predicted trajectory for a single individual. This is why the default
  `start` time was changed for version 0.4 of flexsurv - this was
  previously a vector of the start times observed in the data.

- censtime:

  A right-censoring time, or vector of times matching the rows of
  `newdata`. If `NULL` (the default) then uncensored times to events are
  simulated.

- tidy:

  If `TRUE` then a "tidy" or "long"-format data frame is returned, with
  rows defined by combinations of covariates and simulation replicates.
  The simulation replicate is indicated in the column named `i`.

  If `FALSE`, then a data frame is returned with one row per set of
  covariate values, and different columns for different simulation
  replicates. This is the traditional format for \`simulate\` methods in
  base R.

  In either case, the simulated time and indicator for whether the time
  is an event time (rather than a time of right-censoring) are returned
  in different columns.

- ...:

  Other arguments (not currently used).

## Value

A data frame, with format determined by whether `tidy` was specified.

## Examples

``` r
fit <- flexsurvreg(formula = Surv(futime, fustat) ~ rx, data = ovarian, dist="weibull")
fit2 <- flexsurvspline(formula = Surv(futime, fustat) ~ rx, data = ovarian, k=3)
nd = data.frame(rx=1:2)
simulate(fit, seed=1002, newdata=nd)
#>      time_1 event_1
#> 1  575.5959       1
#> 2 2391.6927       1
simulate(fit, seed=1002, newdata=nd, start=500)
#>      time_1 event_1
#> 1  993.9562       1
#> 2 2750.1742       1
simulate(fit2, nsim=3, seed=1002, newdata=nd)
#>     time_1   time_2    time_3 event_1 event_2 event_3
#> 1 457.8966 3564.307  943.0164       1       1       1
#> 2 419.5862  731.334 9729.0033       1       1       1
simulate(fit2, nsim=3, seed=1002, newdata=nd, start=c(500,1000))
#>     time_1   time_2    time_3 event_1 event_2 event_3
#> 1 2081.447 9364.557  4008.034       1       1       1
#> 2 2872.484 5393.244 24816.827       1       1       1
```
