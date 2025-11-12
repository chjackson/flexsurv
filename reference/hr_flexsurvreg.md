# Hazard ratio as a function of time from a parametric survival model

Hazard ratio as a function of time from a parametric survival model

## Usage

``` r
hr_flexsurvreg(
  x,
  newdata = NULL,
  t = NULL,
  start = 0,
  ci = TRUE,
  B = 1000,
  cl = 0.95,
  na.action = na.pass
)
```

## Arguments

- x:

  Object returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

- newdata:

  A data frame with two rows, each specifying a set of covariate values.
  The hazard ratio is calculated as hazard(z2)/hazard(z1), where z1 is
  the first row of `newdata` and z2 is the second row.

  `newdata` must be supplied unless the model `x` includes just one
  covariate. With one covariate, a default is constructed, which defines
  the hazard ratio between the second and first level of the factor (if
  the covariate is a factor), or between a value of 1 and a value of 0
  (if the covariate is numeric).

- t:

  Times to calculate fitted values for. By default, these are the sorted
  unique observation (including censoring) times in the data - for
  left-truncated datasets these are the "stop" times.

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

- ci:

  Set to `FALSE` to omit confidence intervals.

- B:

  Number of simulations from the normal asymptotic distribution of the
  estimates used to calculate confidence intervals or standard errors.
  Decrease for greater speed at the expense of accuracy, or set `B=0` to
  turn off calculation of CIs and SEs.

- cl:

  Width of symmetric confidence intervals, relative to 1.

- na.action:

  Function determining what should be done with missing values in
  `newdata`. If `na.pass` (the default) then summaries of `NA` are
  produced for missing covariate values. If `na.omit`, then missing
  values are dropped, the behaviour of `summary.flexsurvreg` before
  `flexsurv` version 1.2.

## Value

A data frame with estimate and confidence limits for the hazard ratio,
and one row for each of the times requested in `t`.
