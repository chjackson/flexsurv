# Piecewise exponential distribution for use with flexsurvreg()

Piecewise exponential distribution for use with flexsurvreg()

## Usage

``` r
custom_pexp(cuts, baseline = 1)
```

## Arguments

- cuts:

  Times of the change points (positive, any order). These define
  `K = length(cuts) + 1` periods,
  `[0, cuts[1]], (cuts[1], cuts[2]], ..., (cuts[K-1], Inf)`.

- baseline:

  Index of the period whose hazard is the baseline hazard `rate`. The
  other parameters are hazard ratios named by period, so that the hazard
  in period `k` is `rate * hrk`. For example with K = 3 and
  `baseline = 1` (the default) the parameters are `rate`, `hr2`, `hr3`;
  with `baseline = 2` they are `rate`, `hr1`, `hr3`. `rate` is the
  location parameter, so covariates on it give a proportional hazards
  model.

## Value

A list with components `dist`: the custom distribution list, to pass as
the `dist` argument of `flexsurvreg`, and `dfns`: the list of d, p, h,
H, q, r functions, to pass as the `dfns` argument.

## Examples

``` r
pe <- custom_pexp(c(365, 730))
flexsurvreg(Surv(recyrs*365, censrec) ~ 1, data = bc,
            dist = pe$dist, dfns = pe$dfns)
#> Forming integrated rmst function...
#> Forming integrated mean function...
#> Call:
#> flexsurvreg(formula = Surv(recyrs * 365, censrec) ~ 1, data = bc, 
#>     dist = pe$dist, dfns = pe$dfns)
#> 
#> Estimates: 
#>       est       L95%      U95%      se      
#> rate  2.35e-04  1.81e-04  3.06e-04  3.15e-05
#> hr2   2.40e+00  1.74e+00  3.31e+00  3.94e-01
#> hr3   1.67e+00  1.22e+00  2.29e+00  2.66e-01
#> 
#> N = 686,  Events: 299,  Censored: 387
#> Total time at risk: 771400
#> Log-likelihood = -2632.685, df = 3
#> AIC = 5271.37
#> 
```
