# Generic function to find restricted mean survival time for some distribution

Generic function to find the restricted mean of a distribution, given
the equivalent probability distribution function, using numeric
integration.

## Usage

``` r
rmst_generic(pdist, t, start = 0, matargs = NULL, scalarargs = NULL, ...)
```

## Arguments

- pdist:

  Probability distribution function, for example,
  [`pnorm`](https://rdrr.io/r/stats/Normal.html) for the normal
  distribution, which must be defined in the current workspace. This
  should accept and return vectorised parameters and values. It should
  also return the correct values for the entire real line, for example a
  positive distribution should have `pdist(x)==0` for \\x\<0\\.

- t:

  Vector of times at which rmst is evaluated

- start:

  Optional left-truncation time or times. The returned restricted mean
  survival will be conditioned on survival up to this time.

- matargs:

  Character vector giving the elements of `...` which represent vector
  parameters of the distribution. Empty by default. When vectorised,
  these will become matrices. This is used for the arguments `gamma` and
  `knots` in
  [`psurvspline`](http://chjackson.github.io/flexsurv-dev/reference/Survspline.md).

- scalarargs:

  Character vector naming scalar arguments of the distribution function
  that cannot be vectorised. This is used, for example, for the
  arguments `scale` and `timescale` in
  [`psurvspline`](http://chjackson.github.io/flexsurv-dev/reference/Survspline.md).

- ...:

  The remaining arguments define parameters of the distribution `pdist`.
  These MUST be named explicitly.

## Value

Vector of restricted mean survival times of the distribution at `p`.

## Details

This function is used by default for custom distributions for which an
`rmst` function is not provided.

This assumes a suitably smooth, continuous distribution.

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>

## Examples

``` r
rmst_lnorm(500, start=250, meanlog=7.4225, sdlog = 1.1138)
#> [1] 237.8849
rmst_generic(plnorm, 500, start=250, meanlog=7.4225, sdlog = 1.1138)
#> [1] 237.8849
# must name the arguments
```
