# Bayesian Information Criterion (BIC) for comparison of flexsurvreg models

Bayesian Information Criterion (BIC) for comparison of flexsurvreg
models

## Usage

``` r
# S3 method for class 'flexsurvreg'
BIC(object, cens = TRUE, ...)
```

## Arguments

- object:

  Fitted model returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  (or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)).

- cens:

  Include censored observations in the sample size term (`n`) used in
  the calculation of BIC.

- ...:

  Other arguments (currently unused).

## Value

The BIC of the fitted model. This is minus twice the log likelihood plus
`p*log(n)`, where `p` is the number of parameters and `n` is the sample
size of the data. If `weights` was supplied to `flexsurv`, the sample
size is defined as the sum of the weights.

## Details

There is no "official" definition of what the sample size should be for
the use of BIC in censored survival analysis. BIC is based on an
approximation to Bayesian model comparison using Bayes factors and an
implicit vague prior with the same amount of information as one
observation "unit". Informally, the sample size describes the number of
"units" giving rise to a distinct piece of information (Kass and Raftery
1995). However censored observations provide less information than
observed events, in principle. The default used here is the number of
individuals, for consistency with more familiar kinds of statistical
modelling. However if censoring is heavy, then the number of events may
be a better represent the amount of information. Following these
principles, the best approximation would be expected to be somewhere in
between.

AIC and BIC are intended to measure different things. Briefly, AIC
measures predictive ability, whereas BIC is expected to choose the true
model from a set of models where one of them is the truth. Therefore BIC
chooses simpler models for all but the tiniest sample sizes
(\\log(n)\>2\\, \\n\>7\\). AIC might be preferred in the typical
situation where "all models are wrong but some are useful". AIC also
gives similar results to cross-validation (Stone 1977).

## References

Kass, R. E., & Raftery, A. E. (1995). Bayes factors. Journal of the
American Statistical Association, 90(430), 773-795.

Stone, M. (1977). An asymptotic equivalence of choice of model by
cross‐validation and Akaike's criterion. Journal of the Royal
Statistical Society: Series B (Methodological), 39(1), 44-47.

## See also

[`BIC`](https://rdrr.io/r/stats/AIC.html),
[`AIC`](https://rdrr.io/r/stats/AIC.html),
[`AICC.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/AICc.flexsurvreg.md),
[`nobs.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/nobs.flexsurvreg.md)
