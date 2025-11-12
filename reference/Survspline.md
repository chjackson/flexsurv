# Royston/Parmar spline survival distribution

Probability density, distribution, quantile, random generation, hazard,
cumulative hazard, mean and restricted mean functions for the
Royston/Parmar spline model. These functions have all parameters of the
distribution collected together in a single argument `gamma`. For the
equivalent functions with one argument per parameter, see
[`Survsplinek`](http://chjackson.github.io/flexsurv-dev/reference/Survsplinek.md).

## Usage

``` r
dsurvspline(
  x,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0,
  log = FALSE
)

psurvspline(
  q,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline(
  p,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

rsurvspline(
  n,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0
)

Hsurvspline(
  x,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0
)

hsurvspline(
  x,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0
)

rmst_survspline(
  t,
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0,
  start = 0
)

mean_survspline(
  gamma,
  beta = 0,
  X = 0,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  offset = 0
)
```

## Arguments

- x, q, t:

  Vector of times.

- gamma:

  Parameters describing the baseline spline function, as described in
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).
  This may be supplied as a vector with number of elements equal to the
  length of `knots`, in which case the parameters are common to all
  times. Alternatively a matrix may be supplied, with rows corresponding
  to different times, and columns corresponding to `knots`.

- beta:

  Vector of covariate effects. Not supported and ignored since version
  2.3, and this argument will be removed in 2.4.

- X:

  Matrix of covariate values. Not supported and ignored since version
  2.3, and this argument will be removed in 2.4.

- knots:

  Locations of knots on the axis of log time, supplied in increasing
  order. Unlike in
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  these include the two boundary knots. If there are no additional
  knots, the boundary locations are not used. If there are one or more
  additional knots, the boundary knots should be at or beyond the
  minimum and maximum values of the log times. In
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
  these are exactly at the minimum and maximum values.

  This may in principle be supplied as a matrix, in the same way as for
  `gamma`, but in most applications the knots will be fixed.

- scale:

  `"hazard"`, `"odds"`, or `"normal"`, as described in
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).
  With the default of no knots in addition to the boundaries, this model
  reduces to the Weibull, log-logistic and log-normal respectively. The
  scale must be common to all times.

- timescale:

  `"log"` or `"identity"` as described in
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

- spline:

  `"rp"` to use the natural cubic spline basis described in Royston and
  Parmar. `"splines2ns"` to use the alternative natural cubic spline
  basis from the `splines2` package (Wang and Yan 2021), which may be
  better behaved due to the basis being orthogonal.

- offset:

  An extra constant to add to the linear predictor \\\eta\\. Not
  supported and ignored since version 2.3, and this argument will be
  removed in 2.4.

- log, log.p:

  Return log density or probability.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(X \le x)\\,
  otherwise, \\P(X \> x)\\.

- p:

  Vector of probabilities.

- n:

  Number of random numbers to simulate.

- start:

  Optional left-truncation time or times. The returned restricted mean
  survival will be conditioned on survival up to this time.

## Value

`dsurvspline` gives the density, `psurvspline` gives the distribution
function, `hsurvspline` gives the hazard and `Hsurvspline` gives the
cumulative hazard, as described in
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

`qsurvspline` gives the quantile function, which is computed by crude
numerical inversion (using
[`qgeneric`](http://chjackson.github.io/flexsurv-dev/reference/qgeneric.md)).

`rsurvspline` generates random survival times by using `qsurvspline` on
a sample of uniform random numbers. Due to the numerical root-finding
involved in `qsurvspline`, it is slow compared to typical random number
generation functions.

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

## Examples

``` r
## reduces to the weibull
regscale <- 0.786; cf <- 1.82
a <- 1/regscale; b <- exp(cf)
dweibull(1, shape=a, scale=b)
#> [1] 0.1137858
dsurvspline(1, gamma=c(log(1 / b^a), a)) # should be the same
#> [1] 0.1137858

## reduces to the log-normal
meanlog <- 1.52; sdlog <- 1.11
dlnorm(1, meanlog, sdlog) 
#> [1] 0.1407338
dsurvspline(1, gamma = c(-meanlog/sdlog, 1/sdlog), scale="normal")
#> [1] 0.1407338
# should be the same
```
