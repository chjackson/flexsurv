# Mean time to final state in a mixture multi-state model

Calculate the mean time from the start of the process to a final (or
"absorbing") state in a mixture multi-state model. Models with cycles
are not supported.

## Usage

``` r
meanfinal_fmixmsm(x, newdata = NULL, final = FALSE, B = NULL)
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

## Value

A data frame of mean times to absorption, by covariate values and
pathway (or by final state)
