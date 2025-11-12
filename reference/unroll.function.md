# Convert a function with matrix arguments to a function with vector arguments.

Given a function with matrix arguments, construct an equivalent function
which takes vector arguments defined by the columns of the matrix. The
new function simply uses `cbind` on the vector arguments to make a
matrix, and calls the old one.

## Usage

``` r
unroll.function(mat.fn, ...)
```

## Arguments

- mat.fn:

  A function with any number of arguments, some of which are matrices.

- ...:

  A series of other arguments. Their names define which arguments of
  `mat.fn` are matrices. Their values define a vector of strings to be
  appended to the names of the arguments in the new function. For
  example

  `fn <- unroll.function(oldfn, gamma=1:3, alpha=0:1)`

  will make a new function `fn` with arguments
  `gamma1`,`gamma2`,`gamma3`,`alpha0`,`alpha1`.

  Calling

  `fn(gamma1=a,gamma2=b,gamma3=c,alpha0=d,alpha1=e)`

  should give the same answer as

  `oldfn(gamma=cbind(a,b,c),alpha=cbind(d,e))`

## Value

The new function, with vector arguments.

## Usage in flexsurv

This is used by
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
to allow spline models, which have an arbitrary number of parameters, to
be fitted using
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

The “custom distributions” facility of
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
expects the user-supplied probability density and distribution functions
to have one explicitly named argument for each scalar parameter, and
given R vectorisation, each of those arguments could be supplied as a
vector of alternative parameter values.

However, spline models have a varying number of scalar parameters,
determined by the number of knots in the spline.
[`dsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/Survspline.md)
and
[`psurvspline`](http://chjackson.github.io/flexsurv-dev/reference/Survspline.md)
have an argument called `gamma`. This can be supplied as a matrix, with
number of columns `n` determined by the number of knots (plus 2), and
rows referring to alternative parameter values. The following statements
are used in the source of `flexsurvspline`:

     dfn <-
    unroll.function(dsurvspline, gamma=0:(nk-1)) pfn <-
    unroll.function(psurvspline, gamma=0:(nk-1)) 

to convert these into functions with arguments `gamma0`,
`gamma1`,...,`gamman`, corresponding to the columns of `gamma`, where
`n = nk-1`, and with other arguments in the same format.

## See also

[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>

## Examples

``` r
fn <- unroll.function(ncol, x=1:3)
fn(1:3, 1:3, 1:3) # equivalent to...
#> [1] 3
ncol(cbind(1:3,1:3,1:3))
#> [1] 3
```
