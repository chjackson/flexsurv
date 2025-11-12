# Transition-specific parameters in a flexible parametric multi-state model

List of maximum likelihood estimates of transition-specific parameters
in a flexible parametric multi-state model, at given covariate values.

## Usage

``` r
pars.fmsm(x, trans, newdata = NULL, tvar = "trans")
```

## Arguments

- x:

  A multi-state model fitted with
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md)
  for the required form of the model and the data.

  `x` can also be a list of
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  models, with one component for each permitted transition in the
  multi-state model, as illustrated in
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- trans:

  Matrix indicating allowed transitions. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- newdata:

  A data frame specifying the values of covariates in the fitted model,
  other than the transition number. See
  [`msfit.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/msfit.flexsurvreg.md).

- tvar:

  Variable in the data representing the transition type. Not required if
  `x` is a list of models.

## Value

A list with one component for each permitted transition. Each component
has one element for each parameter of the parametric distribution that
generates the corresponding event in the multi-state model.

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>.
