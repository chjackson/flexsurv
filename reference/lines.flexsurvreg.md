# Add fitted flexible survival curves to a plot

Add fitted survival (or hazard or cumulative hazard) curves from a
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
model fit to an existing plot.

## Usage

``` r
# S3 method for class 'flexsurvreg'
lines(
  x,
  newdata = NULL,
  X = NULL,
  type = "survival",
  t = NULL,
  est = TRUE,
  ci = NULL,
  B = 1000,
  cl = 0.95,
  col = "red",
  lty = 1,
  lwd = 2,
  col.ci = NULL,
  lty.ci = 2,
  lwd.ci = 1,
  ...
)
```

## Arguments

- x:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
  representing a fitted survival model object.

- newdata:

  Covariate values to produce fitted curves for, as a data frame, as
  described in
  [`plot.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/plot.flexsurvreg.md).

- X:

  Covariate values to produce fitted curves for, as a matrix, as
  described in
  [`plot.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/plot.flexsurvreg.md).

- type:

  `"survival"` for survival, `"cumhaz"` for cumulative hazard, or
  `"hazard"` for hazard, as in
  [`plot.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/plot.flexsurvreg.md).

- t:

  Vector of times to plot fitted values for.

- est:

  Plot fitted curves (`TRUE` or `FALSE`.)

- ci:

  Plot confidence intervals for fitted curves.

- B:

  Number of simulations controlling accuracy of confidence intervals, as
  used in
  [`summary`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md).

- cl:

  Width of confidence intervals, by default 0.95 for 95% intervals.

- col:

  Colour of the fitted curve(s).

- lty:

  Line type of the fitted curve(s).

- lwd:

  Line width of the fitted curve(s).

- col.ci:

  Colour of the confidence limits, defaulting to the same as for the
  fitted curve.

- lty.ci:

  Line type of the confidence limits.

- lwd.ci:

  Line width of the confidence limits, defaulting to the same as for the
  fitted curve.

- ...:

  Other arguments to be passed to the generic
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) and
  [`lines`](https://rdrr.io/r/graphics/lines.html) functions.

## Details

Equivalent to
[`plot.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/plot.flexsurvreg.md)`(...,add=TRUE)`.

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
