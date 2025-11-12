# Royston/Parmar spline survival distribution functions with one argument per parameter

Probability density, distribution, quantile, random generation, hazard,
cumulative hazard, mean and restricted mean functions for the
Royston/Parmar spline model, with one argument per parameter. For the
equivalent functions with all parameters collected together in a single
argument, see
[`Survspline`](http://chjackson.github.io/flexsurv-dev/reference/Survspline.md).

## Usage

``` r
mean_survspline0(
  gamma0,
  gamma1,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline1(
  gamma0,
  gamma1,
  gamma2,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline2(
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline3(
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline4(
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline5(
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline6(
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

mean_survspline7(
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rmst_survspline0(
  t,
  gamma0,
  gamma1,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline1(
  t,
  gamma0,
  gamma1,
  gamma2,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline2(
  t,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline3(
  t,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline4(
  t,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline5(
  t,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline6(
  t,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

rmst_survspline7(
  t,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots = c(-10, 10),
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  start = 0
)

dsurvspline0(
  x,
  gamma0,
  gamma1,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline1(
  x,
  gamma0,
  gamma1,
  gamma2,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline2(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline3(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline4(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline5(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline6(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

dsurvspline7(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  log = FALSE
)

psurvspline0(
  q,
  gamma0,
  gamma1,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline1(
  q,
  gamma0,
  gamma1,
  gamma2,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline2(
  q,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline3(
  q,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline4(
  q,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline5(
  q,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline6(
  q,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

psurvspline7(
  q,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline0(
  p,
  gamma0,
  gamma1,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline1(
  p,
  gamma0,
  gamma1,
  gamma2,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline2(
  p,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline3(
  p,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline4(
  p,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline5(
  p,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline6(
  p,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

qsurvspline7(
  p,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp",
  lower.tail = TRUE,
  log.p = FALSE
)

rsurvspline0(
  n,
  gamma0,
  gamma1,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline1(
  n,
  gamma0,
  gamma1,
  gamma2,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline2(
  n,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline3(
  n,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline4(
  n,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline5(
  n,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline6(
  n,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

rsurvspline7(
  n,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline0(
  x,
  gamma0,
  gamma1,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline1(
  x,
  gamma0,
  gamma1,
  gamma2,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline2(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline3(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline4(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline5(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline6(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

hsurvspline7(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline0(
  x,
  gamma0,
  gamma1,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline1(
  x,
  gamma0,
  gamma1,
  gamma2,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline2(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline3(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline4(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline5(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline6(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)

Hsurvspline7(
  x,
  gamma0,
  gamma1,
  gamma2,
  gamma3,
  gamma4,
  gamma5,
  gamma6,
  gamma7,
  gamma8,
  knots,
  scale = "hazard",
  timescale = "log",
  spline = "rp"
)
```

## Arguments

- gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7,
  gamma8:

  Parameters describing the baseline spline function, as described in
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

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

- start:

  Optional left-truncation time or times. The returned restricted mean
  survival will be conditioned on survival up to this time.

- x, q, t:

  Vector of times.

- log, log.p:

  Return log density or probability.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(X \le x)\\,
  otherwise, \\P(X \> x)\\.

- p:

  Vector of probabilities.

- n:

  Number of random numbers to simulate.

## Details

These functions go up to 7 spline knots, or 9 parameters.

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
