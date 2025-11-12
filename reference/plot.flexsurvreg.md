# Plots of fitted flexible survival models

Plot fitted survival, cumulative hazard or hazard from a parametric
model against nonparametric estimates to diagnose goodness-of-fit.
Alternatively plot a user-defined function of the model parameters
against time.

## Usage

``` r
# S3 method for class 'flexsurvreg'
plot(
  x,
  newdata = NULL,
  X = NULL,
  type = "survival",
  fn = NULL,
  t = NULL,
  start = 0,
  est = TRUE,
  ci = NULL,
  B = 1000,
  cl = 0.95,
  col.obs = "black",
  lty.obs = 1,
  lwd.obs = 1,
  col = "red",
  lty = 1,
  lwd = 2,
  col.ci = NULL,
  lty.ci = 2,
  lwd.ci = 1,
  ylim = NULL,
  add = FALSE,
  ...
)
```

## Arguments

- x:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- newdata:

  Data frame containing covariate values to produce fitted values for.
  See
  [`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md).

  If there are only factor covariates in the model, then Kaplan-Meier
  (or nonparametric hazard...) curves are plotted for all distinct
  groups, and by default, fitted curves are also plotted for these
  groups. To plot Kaplan-Meier and fitted curves for only a subset of
  groups, use `plot(survfit())` followed by
  [`lines.flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/lines.flexsurvreg.md).

  If there are any continuous covariates, then a single population
  Kaplan-Meier curve is drawn. By default, a single fitted curve is
  drawn with the covariates set to their mean values in the data - for
  categorical covariates, the means of the 0/1 indicator variables are
  taken.

- X:

  Alternative way to supply covariate values, as a model matrix. See
  [`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md).
  `newdata` is an easier way.

- type:

  `"survival"` for survival, to be plotted against Kaplan-Meier
  estimates from
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html).

  `"cumhaz"` for cumulative hazard, plotted against transformed
  Kaplan-Meier estimates from
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html).

  `"hazard"` for hazard, to be plotted against smooth nonparametric
  estimates from [`muhaz`](https://rdrr.io/pkg/muhaz/man/muhaz.html).
  The nonparametric estimates tend to be unstable, and these plots are
  intended just to roughly indicate the shape of the hazards through
  time. The `min.time` and `max.time` options to
  [`muhaz`](https://rdrr.io/pkg/muhaz/man/muhaz.html) may sometimes need
  to be passed as arguments to `plot.flexsurvreg` to avoid an error
  here.

  Ignored if `"fn"` is specified.

- fn:

  Custom function of the parameters to summarise against time. The first
  two arguments of the function must be `t` representing time, and
  `start` representing left-truncation points, and any remaining
  arguments must be parameters of the distribution. It should return a
  vector of the same length as `t`.

- t:

  Vector of times to plot fitted values for, see
  [`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md).

- start:

  Left-truncation points, see
  [`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md).

- est:

  Plot fitted curves (`TRUE` or `FALSE`.)

- ci:

  Plot confidence intervals for fitted curves. By default, this is
  `TRUE` if one observed/fitted curve is plotted, and `FALSE` if
  multiple curves are plotted.

- B:

  Number of simulations controlling accuracy of confidence intervals, as
  used in
  [`summary`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md).
  Decrease for greater speed at the expense of accuracy, or set `B=0` to
  turn off calculation of CIs.

- cl:

  Width of confidence intervals, by default 0.95 for 95% intervals.

- col.obs:

  Colour of the nonparametric curve.

- lty.obs:

  Line type of the nonparametric curve.

- lwd.obs:

  Line width of the nonparametric curve.

- col:

  Colour of the fitted parametric curve(s).

- lty:

  Line type of the fitted parametric curve(s).

- lwd:

  Line width of the fitted parametric curve(s).

- col.ci:

  Colour of the fitted confidence limits, defaulting to the same as for
  the fitted curve.

- lty.ci:

  Line type of the fitted confidence limits.

- lwd.ci:

  Line width of the fitted confidence limits.

- ylim:

  y-axis limits: vector of two elements.

- add:

  If `TRUE`, add lines to an existing plot, otherwise new axes are
  drawn.

- ...:

  Other options to be passed to
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html)
  or [`muhaz`](https://rdrr.io/pkg/muhaz/man/muhaz.html), for example,
  to control the smoothness of the nonparametric hazard estimates. The
  `min.time` and `max.time` options to
  [`muhaz`](https://rdrr.io/pkg/muhaz/man/muhaz.html) may sometimes need
  to be changed from the defaults.

## Note

Some standard plot arguments such as `"xlim","xlab"` may not work. This
function was designed as a quick check of model fit. Users wanting
publication-quality graphs are advised to set up an empty plot with the
desired axes first (e.g. with `plot(...,type="n",...)`), then use
suitable [`lines`](https://rdrr.io/r/graphics/lines.html) functions to
add lines.

If case weights were used to fit the model, these are used when
producing nonparametric estimates of survival and cumulative hazard, but
not for the hazard estimates.

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
