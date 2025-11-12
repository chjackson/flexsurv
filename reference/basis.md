# Natural cubic spline basis

Compute a basis for a natural cubic spline, by default using the
parameterisation described by Royston and Parmar (2002). Used for
flexible parametric survival models.

## Usage

``` r
basis(knots, x, spline = "rp")
```

## Arguments

- knots:

  Vector of knot locations in increasing order, including the boundary
  knots at the beginning and end.

- x:

  Vector of ordinates to compute the basis for.

- spline:

  `"rp"` to use the natural cubic spline basis described in Royston and
  Parmar. `"splines2ns"` to use the alternative natural cubic spline
  basis from the `splines2` package (Wang and Yan 2021), which may be
  better behaved due to the basis being orthogonal.

## Value

A matrix with one row for each ordinate and one column for each knot.

`basis` returns the basis, and `dbasis` returns its derivative with
respect to `x`.

`fss` and `dfss` are the same, but with the order of the arguments
swapped around for consistency with similar functions in other R
packages.

## Details

The exact formula for the basis is given in
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## References

Royston, P. and Parmar, M. (2002). Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of
treatment effects. Statistics in Medicine 21(1):2175-2197.

Wang W, Yan J (2021). Shape-Restricted Regression Splines with R Package
splines2. Journal of Data Science, 19(3), 498-517.

## See also

[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
