# Summarise quantities of interest from fitted flexsurvrtrunc models

This function extracts quantities of interest from the untruncated
version of a model with individual-specific right truncation points
fitted by
[`flexsurvrtrunc`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvrtrunc.md).
Note that covariates are currently not supported by
[`flexsurvrtrunc`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvrtrunc.md).

## Usage

``` r
# S3 method for class 'flexsurvrtrunc'
summary(
  object,
  type = "survival",
  fn = NULL,
  t = NULL,
  quantiles = 0.5,
  ci = TRUE,
  se = FALSE,
  B = 1000,
  cl = 0.95,
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

- type:

  `"survival"` for survival probabilities.

  `"cumhaz"` for cumulative hazards.

  `"hazard"` for hazards.

  `"rmst"` for restricted mean survival.

  `"mean"` for mean survival.

  `"median"` for median survival (alternative to `type="quantile"` with
  `quantiles=0.5`).

  `"quantile"` for quantiles of the survival time distribution.

  Ignored if `"fn"` is specified.

- fn:

  Custom function of the parameters to summarise against time. This has
  optional first argument `t` representing time, and any remaining
  arguments must be parameters of the distribution. It should return a
  vector of the same length as `t`.

- t:

  Times to calculate fitted values for. By default, these are the sorted
  unique observation (including censoring) times in the data - for
  left-truncated datasets these are the "stop" times.

- quantiles:

  If `type="quantile"`, this specifies the quantiles of the survival
  time distribution to return estimates for.

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

- ...:

  Further arguments passed to or from other methods. Currently unused.
