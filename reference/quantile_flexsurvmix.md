# Quantiles of time-to-event distributions in a flexsurvmix model

This returns the quantiles of each event-specific parametric
time-to-event distribution in the mixture model, which describes the
time to the event conditionally on that event being the one that
happens.

## Usage

``` r
quantile_flexsurvmix(x, newdata = NULL, B = NULL, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- x:

  Fitted model object returned from
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).

- newdata:

  Data frame or list of covariate values. If omitted for a model with
  covariates, a default is used, defined by all combinations of factors
  if the only covariates in the model are factors, or all covariate
  values of zero if there are any non-factor covariates in the model.

- B:

  Number of simulations to use to compute 95% confidence intervals,
  based on the asymptotic multivariate normal distribution of the basic
  parameter estimates. If `B=NULL` then intervals are not computed.

- probs:

  Vector of alternative quantiles, by default `c(0.025, 0.95, 0.975)`
  giving the median and a 95% interval.
