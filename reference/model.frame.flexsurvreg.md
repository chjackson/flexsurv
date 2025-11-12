# Extract original data from `flexsurvreg` objects.

Extract the data from a model fitted with `flexsurvreg`.

## Usage

``` r
# S3 method for class 'flexsurvreg'
model.frame(formula, ...)

# S3 method for class 'flexsurvreg'
model.matrix(object, par = NULL, ...)
```

## Arguments

- formula:

  A fitted model object, as returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- ...:

  Further arguments (not used).

- object:

  A fitted model object, as returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- par:

  String naming the parameter whose linear model matrix is desired.

  The default value of `par=NULL` returns a matrix consisting of the
  model matrices for all models in the object `cbind`ed together, with
  the intercepts excluded. This is not really a “model matrix” in the
  usual sense, however, the columns directly correspond to the covariate
  coefficients in the matrix of estimates from the fitted model.

## Value

`model.frame` returns a data frame with all the original variables used
for the model fit.

`model.matrix` returns a design matrix for a part of the model that
includes covariates. The required part is indicated by the `"par"`
argument (see above).

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`model.frame`](https://rdrr.io/r/stats/model.frame.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
