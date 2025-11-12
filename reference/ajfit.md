# Aalen-Johansen nonparametric estimates comparable to a fitted flexsurvmix model

Given a fitted flexsurvmix model, return the Aalen-Johansen estimates of
the probability of occupying each state at a series of times covering
the observed data. State 1 represents not having experienced any of the
competing events, while state 2 and any further states correspond to
having experienced each of the competing events respectively. These
estimates can be compared with the fitted probabilities returned by
[`p_flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/p_flexsurvmix.md)
to check the fit of a `flexsurvmix` model.

## Usage

``` r
ajfit(x, newdata = NULL, tidy = TRUE)
```

## Arguments

- x:

  Fitted model returned by
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).

- newdata:

  Data frame of alternative covariate values to check fit for. Only
  factor covariates are supported.

- tidy:

  If `TRUE` then a single tidy data frame is returned. Otherwise the
  function returns the object returned by `survfit`, or a list of these
  objects if we are computing subset-specific estimates.

## Details

This is only supported for models with no covariates or models
containing only factor covariates.

For models with factor covariates, the Aalen-Johansen estimates are
computed for the subsets of the data defined in `newdata`. If `newdata`
is not supplied, then this function returns state occupancy
probabilities for all possible combinations of the factor levels.

The Aalen-Johansen estimates are computed using
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) from the
`survival` package (Therneau 2020).

## References

Therneau T (2020). \_A Package for Survival Analysis in R\_. R package
version 3.2-3, \<URL: https://CRAN.R-project.org/package=survival\>.
