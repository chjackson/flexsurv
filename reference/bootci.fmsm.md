# Bootstrap confidence intervals for flexsurv output functions

Calculate a confidence interval for a model output by repeatedly
replacing the parameters in a fitted model object with a draw from the
multivariate normal distribution of the maximum likelihood estimates,
then recalculating the output function.

## Usage

``` r
bootci.fmsm(
  x,
  B,
  fn,
  cl = 0.95,
  attrs = NULL,
  cores = NULL,
  sample = FALSE,
  ...
)
```

## Arguments

- x:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object. Or a list of such
  objects, defining a multi-state model.

- B:

  Number of parameter draws to use

- fn:

  Function to bootstrap the results of. It must have an argument named
  `x` giving a fitted flexsurv model object. This may return a value
  with any format, e.g. list, matrix or vector, as long as it can be
  converted to a numeric vector with `unlist`. See the example below.

- cl:

  Width of symmetric confidence interval, by default 0.95

- attrs:

  Any attributes of the value returned from `fn` which we want
  confidence intervals for. These will be unlisted, if possible, and
  appended to the result vector.

- cores:

  Number of cores to use for parallel processing.

- sample:

  If `TRUE` then the bootstrap sample itself is returned. If `FALSE`
  then the quantiles of the sample are returned giving a confidence
  interval.

- ...:

  Additional arguments to pass to `fn`.

## Value

A matrix with two rows, giving the upper and lower confidence limits
respectively. Each row is a vector of the same length as the unlisted
result of the function corresponding to `fncall`.

## Examples

``` r
## How to use bootci.fmsm

## Write a function with one argument called x giving a fitted model,
## and returning some results of the model.  The results may be in any form.   

tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")

summfn <- function(x, t){
 resp <-  flexsurv::pmatrix.fs(x, trans=tmat, t=t)
 rest <- flexsurv::totlos.fs(x, trans=tmat, t=t)
 list(resp, rest)
}

## Use bootci.fmsm to obtain the confidence interval
## The matrix columns are in the order of the unlisted results of the original
## summfn.  You will have to rearrange them into the format that you want.
## If summfn has any extra arguments, in this case \code{t}, make sure they are
## passed through via the ... argument to bootci.fmsm

bootci.fmsm(bexp, B=3, fn=summfn, t=10)
#>             [,1] [,2] [,3]      [,4]       [,5] [,6]      [,7]      [,8] [,9]
#> 2.5%  0.06379407    0    0 0.1232301 0.05183498    0 0.7364390 0.9118225    1
#> 97.5% 0.09767351    0    0 0.1698328 0.08817746    0 0.8127567 0.9481650    1
#>          [,10] [,11] [,12]    [,13]    [,14] [,15]    [,16]    [,17] [,18]
#> 2.5%  3.399433     0     0 1.955466 3.203392     0 3.926742 6.247245    10
#> 97.5% 3.879065     0     0 2.249966 3.752755     0 4.554488 6.796608    10
bootci.fmsm(bexp, B=3, fn=summfn, t=5)
#>            [,1] [,2] [,3]      [,4]      [,5] [,6]      [,7]      [,8] [,9]
#> 2.5%  0.2949335    0    0 0.2551559 0.2632568    0 0.4145488 0.7130623    1
#> 97.5% 0.3158727    0    0 0.2714704 0.2869377    0 0.4459622 0.7367432    1
#>          [,10] [,11] [,12]    [,13]    [,14] [,15]     [,16]    [,17] [,18]
#> 2.5%  2.887233     0     0 1.020261 2.760045     0 0.9831934 2.144306     5
#> 97.5% 2.968191     0     0 1.072831 2.855694     0 1.0761249 2.239955     5
```
