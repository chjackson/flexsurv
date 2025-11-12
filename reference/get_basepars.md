# Evaluate baseline time-to-event distribution parameters given covariate values in a flexsurvmix model

Evaluate baseline time-to-event distribution parameters given covariate
values in a flexsurvmix model

## Usage

``` r
get_basepars(x, newdata, event)
```

## Arguments

- x:

  Fitted model object

- newdata:

  Data frame of alternative covariate values

- event:

  Event

## Value

A list with one component per parameter, that gives values of that
parameter for the different covariate values
