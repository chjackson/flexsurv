# Quantiles of the distribution of the time until reaching a final state in a mixture multi-state model

Calculate the quantiles of the time from the start of the process to
each possible final (or "absorbing") state in a mixture multi-state
model. Models with cycles are not supported.

## Usage

``` r
qfinal_fmixmsm(
  x,
  newdata = NULL,
  final = FALSE,
  B = NULL,
  n = 10000,
  probs = c(0.025, 0.5, 0.975)
)
```

## Arguments

- x:

  Object returned by
  [`fmixmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmixmsm.md),
  representing a multi-state model built from piecing together mixture
  models fitted by
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).

- newdata:

  Data frame or list of covariate values. If omitted for a model with
  covariates, a default is used, defined by all combinations of factors
  if the only covariates in the model are factors, or all covariate
  values of zero if there are any non-factor covariates in the model.

- final:

  If `TRUE` then the mean time to the final state is calculated for each
  final state, by taking a weighted average of the mean time to travel
  each pathway ending in that final state, weighted by the probability
  of the pathway. If `FALSE` (the default) then a separate mean is
  calculated for each pathway.

- B:

  Number of simulations to use to compute 95% confidence intervals,
  based on the asymptotic multivariate normal distribution of the basic
  parameter estimates. If `B=NULL` then intervals are not computed.

- n:

  Number of individual-level simulations to use to characterise the
  time-to-event distributions

- probs:

  Quantiles to calculate, by default, `c(0.025, 0.5, 0.975)`

## Value

Data frame of quantiles of the time to final state by pathway and
covariate value, or by final state and covariate value.
