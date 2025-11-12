# Second-order Akaike information criterion

Second-order or "corrected" Akaike information criterion, often known as
AICc (see, e.g. Burnham and Anderson 2002). This is defined as -2
log-likelihood + `(2*p*n)/(n - p -1)`, whereas the standard AIC is
defined as -2 log-likelihood + `2*p`, where `p` is the number of
parameters and `n` is the sample size. The correction is intended to
adjust AIC for small-sample bias, hence it only makes a difference to
the result for small `n`.

## Usage

``` r
# S3 method for class 'flexsurvreg'
AICc(object, cens = TRUE, ...)

# S3 method for class 'flexsurvreg'
AICC(object, cens = TRUE, ...)
```

## Arguments

- object:

  Fitted model returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  (or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)).

- cens:

  Include censored observations in the sample size term (`n`) used in
  this calculation. See
  [`BIC.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/BIC.flexsurvreg.md)
  for a discussion of the issues with defining the sample size.

- ...:

  Other arguments (currently unused).

## Value

The second-order AIC of the fitted model.

## Details

This can be spelt either as `AICC` or `AICc`.

## References

Burnham, K. P., Anderson, D. R. (2002) Model Selection and Multimodel
Inference: a practical information-theoretic approach. Second edition.
Springer: New York.

## See also

[`BIC`](https://rdrr.io/r/stats/AIC.html),
[`AIC`](https://rdrr.io/r/stats/AIC.html),
[`BIC.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/BIC.flexsurvreg.md),
[`nobs.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/nobs.flexsurvreg.md)
