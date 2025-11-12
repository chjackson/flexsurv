# Extract model coefficients from fitted flexible survival models

Extract model coefficients from fitted flexible survival models. This
presents all parameter estimates, transformed to the real line if
necessary. For example, shape or scale parameters, which are constrained
to be positive, are returned on the log scale.

## Usage

``` r
# S3 method for class 'flexsurvreg'
coef(object, ...)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- ...:

  Further arguments passed to or from other methods. Currently unused.

## Value

This returns the `mod$res.t[,"est"]` component of the fitted model
object `mod`. See
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
for full documentation of all components.

## Details

This matches the behaviour of `coef.default` for standard R model
families such as [`glm`](https://rdrr.io/r/stats/glm.html), where
intercepts in regression models are presented on the same scale as the
covariate effects. Note that any parameter in a distribution fitted by
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
or
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
may be an intercept in a regression model.

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
