# Summaries of fitted flexible survival models

Return fitted survival, cumulative hazard or hazard at a series of times
from a fitted
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
or
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
model.

## Usage

``` r
# S3 method for class 'flexsurvreg'
summary(
  object,
  newdata = NULL,
  X = NULL,
  type = "survival",
  fn = NULL,
  t = NULL,
  quantiles = 0.5,
  start = 0,
  cross = TRUE,
  ci = TRUE,
  se = FALSE,
  B = 1000,
  cl = 0.95,
  tidy = FALSE,
  na.action = na.pass,
  ...
)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- newdata:

  Data frame containing covariate values to produce fitted values for.
  Or a list that can be coerced to such a data frame. There must be a
  column for every covariate in the model formula, and one row for every
  combination of covariates the fitted values are wanted for. These are
  in the same format as the original data, with factors as a single
  variable, not 0/1 contrasts.

  If this is omitted, if there are any continuous covariates, then a
  single summary is provided with all covariates set to their mean
  values in the data - for categorical covariates, the means of the 0/1
  indicator variables are taken. If there are only factor covariates in
  the model, then all distinct groups are used by default.

- X:

  Alternative way of defining covariate values to produce fitted values
  for. Since version 0.4, `newdata` is an easier way that doesn't
  require the user to create factor contrasts, but `X` has been kept for
  backwards compatibility.

  Columns of `X` represent different covariates, and rows represent
  multiple combinations of covariate values. For example
  `matrix(c(1,2),nrow=2)` if there is only one covariate in the model,
  and we want survival for covariate values of 1 and 2. A vector can
  also be supplied if just one combination of covariates is needed.

  For “factor” (categorical) covariates, the values of the contrasts
  representing factor levels (as returned by the
  [`contrasts`](https://rdrr.io/r/stats/contrasts.html) function) should
  be used. For example, for a covariate `agegroup` specified as an
  unordered factor with levels `20-29, 30-39, 40-49, 50-59`, and
  baseline level `20-29`, there are three contrasts. To return summaries
  for groups `20-29` and `40-49`, supply
  `X = rbind(c(0,0,0), c(0,1,0))`, since all contrasts are zero for the
  baseline level, and the second contrast is “turned on” for the third
  level `40-49`.

- type:

  `"survival"` for survival probabilities.

  `"cumhaz"` for cumulative hazards.

  `"hazard"` for hazards.

  `"rmst"` for restricted mean survival.

  `"mean"` for mean survival.

  `"median"` for median survival (alternative to `type="quantile"` with
  `quantiles=0.5`).

  `"quantile"` for quantiles of the survival time distribution.

  `"link"` for the fitted value of the location parameter (i.e. the
  "linear predictor" but on the natural scale of the parameter, not on
  the log scale)

  Ignored if `"fn"` is specified.

- fn:

  Custom function of the parameters to summarise against time. This has
  optional first two arguments `t` representing time, and `start`
  representing left-truncation points, and any remaining arguments must
  be parameters of the distribution. It should be vectorised, and return
  a vector corresponding to the vectors given by `t`, `start` and the
  parameter vectors.

- t:

  Times to calculate fitted values for. By default, these are the sorted
  unique observation (including censoring) times in the data - for
  left-truncated datasets these are the "stop" times.

- quantiles:

  If `type="quantile"`, this specifies the quantiles of the survival
  time distribution to return estimates for.

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

- cross:

  If `TRUE` (the default) then summaries are calculated for all
  combinations of times specified in `t` and covariate vectors specifed
  in `newdata`.

  If `FALSE`, then the times `t` should be of length equal to the number
  of rows in `newdata`, and one summary is produced for each row of
  `newdata` paired with the corresponding element of `t`. This is used,
  e.g. when determining Cox-Snell residuals.

- ci:

  Set to `FALSE` to omit confidence intervals.

- se:

  Set to `TRUE` to include standard errors.

- B:

  Number of simulations from the normal asymptotic distribution of the
  estimates used to calculate confidence intervals or standard errors.
  Decrease for greater speed at the expense of accuracy, or set `B=0` to
  turn off calculation of CIs and SEs.

- cl:

  Width of symmetric confidence intervals, relative to 1.

- tidy:

  If `TRUE`, then the results are returned as a tidy data frame instead
  of a list. This can help with using the ggplot2 package to compare
  summaries for different covariate values.

- na.action:

  Function determining what should be done with missing values in
  `newdata`. If `na.pass` (the default) then summaries of `NA` are
  produced for missing covariate values. If `na.omit`, then missing
  values are dropped, the behaviour of `summary.flexsurvreg` before
  `flexsurv` version 1.2.

- ...:

  Further arguments passed to or from other methods. Currently unused.

## Value

If `tidy=FALSE`, a list with one component for each unique covariate
value (if there are only categorical covariates) or one component (if
there are no covariates or any continuous covariates). Each of these
components is a matrix with one row for each time in `t`, giving the
estimated survival (or cumulative hazard, or hazard) and 95% confidence
limits. These list components are named with the covariate names and
values which define them.

If `tidy=TRUE`, a data frame is returned instead. This is formed by
stacking the above list components, with additional columns to identify
the covariate values that each block corresponds to.

If there are multiple summaries, an additional list component named `X`
contains a matrix with the exact values of contrasts (dummy covariates)
defining each summary.

The
[`plot.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/plot.flexsurvreg.md)
function can be used to quickly plot these model-based summaries against
empirical summaries such as Kaplan-Meier curves, to diagnose model fit.

Confidence intervals are obtained by sampling randomly from the
asymptotic normal distribution of the maximum likelihood estimates and
then taking quantiles (see, e.g. Mandel (2013)).

## Details

Time-dependent covariates are not currently supported. The covariate
values are assumed to be constant through time for each fitted curve.

## References

Mandel, M. (2013). "Simulation based confidence intervals for functions
with complicated derivatives." The American Statistician (in press).

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
