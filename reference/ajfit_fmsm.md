# Check the fit of Markov flexible parametric multi-state models against nonparametric estimates.

Computes both parametric and comparable Aalen-Johansen nonparametric
estimates from a flexible parametric multi-state model, and returns them
together in a tidy data frame. Only models with no covariates, or only
factor covariates, are supported. If there are factor covariates, then
the nonparametric estimates are computed for subgroups defined by
combinations of the covariates. Another restriction of this function is
that all transitions must have the same covariates on them.

## Usage

``` r
ajfit_fmsm(x, maxt = NULL, newdata = NULL)
```

## Arguments

- x:

  Object returned by
  [`fmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmsm.md)
  representing a flexible parametric multi-state model. This must be
  Markov, rather than semi-Markov, and no check is performed for this.
  Note that all "competing risks" style models, with just one source
  state and multiple destination states, are Markov, so those are fine
  here.

- maxt:

  Maximum time to compute parametric estimates to.

- newdata:

  Data frame defining the subgroups to consider. This must have a column
  for each covariate in the model. If omitted, then all potential
  subgroups defined by combinations of factor covariates are included.
  Continuous covariates are not supported.

## Value

Tidy data frame containing both Aalen-Johansen and parametric estimates
of state occupancy over time, and by subgroup if subgroups are included.
