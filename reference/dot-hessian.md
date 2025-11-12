# Numerical evaluation of the hessian of a function using numDeriv::hessian

We perform a quick check about the expected runtime and adjust the
precision accordingly.

## Usage

``` r
.hessian(f, x, seconds.warning = 60, default.r = 6, min.r = 2, ...)
```

## Arguments

- f:

  function to compute Hessian for

- x:

  location to evaluate Hessian at

- seconds.warning:

  time threshold in seconds to trigger message and reduce the number of
  iterations for Richardson extrapolation of numDeriv::hessian

- default.r:

  default number of iterations (high-precision recommendation of
  numDeriv)

- min.r:

  minial number of iteration, must be at least 2,

- ...:

  further arguments passed to method.args of numDeriv::hessian
