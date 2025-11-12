# Marginal survival and hazards of fitted flexsurvreg models

Returns a tidy data.frame of marginal survival probabilities, or
hazards, restricted mean survival, or quantiles of the marginal survival
function at user-defined time points and covariate patterns.
Standardization is performed over any undefined covariates in the model.
The user provides the data to standardize over. Contrasts can be
calculated resulting in estimates of the average treatment effect or the
average treatment effect in the treated if a treated subset of the data
are supplied.

## Usage

``` r
standsurv(
  object,
  newdata = NULL,
  at = list(list()),
  atreference = 1,
  type = "survival",
  t = NULL,
  ci = FALSE,
  se = FALSE,
  boot = FALSE,
  B = NULL,
  cl = 0.95,
  trans = "log",
  contrast = NULL,
  trans.contrast = NULL,
  seed = NULL,
  rmap,
  ratetable,
  scale.ratetable = 365.25,
  n.gauss.quad = 100,
  quantiles = 0.5,
  interval = c(1e-08, 500)
)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- newdata:

  Data frame containing covariate values to produce marginal values for.
  If not specified then the fitted model data.frame is used. There must
  be a column for every covariate in the model formula for which the
  user wishes to standardize over. These are in the same format as the
  original data, with factors as a single variable, not 0/1 contrasts.
  Any covariates that are to be fixed should be specified in `at`. There
  should be one row for every combination of covariates in which to
  standardize over. If newdata contains a variable named '(weights)'
  then a weighted mean will be used to create the standardized
  estimates. This is the default behaviour if the fitted model contains
  case weights, which are stored in the fitted model data.frame.

- at:

  A list of scenarios in which specific covariates are fixed to certain
  values. Each element of `at` must itself be a list. For example, for a
  covariate `group` with levels "Good", "Medium" and "Poor", the
  standardized survival plots for each group averaging over all other
  covariates is specified using
  `at=list(list(group="Good"), list(group="Medium"), list(group="Poor"))`.

- atreference:

  The reference scenario for making contrasts. Default is 1 (i.e. the
  first element of `at`).

- type:

  `"survival"` for marginal survival probabilities. In a relative
  survival framework this returns the marginal all-cause survival (see
  details).

  `"hazard"` for the hazard of the marginal survival probability. In a
  relative survival framework this returns the marginal all-cause hazard
  (see details).

  `"rmst"` for standardized restricted mean survival.

  `"relsurvival"` for marginal relative survival (can only be specified
  if a relative survival model has been fitted in flexsurv).

  `"excesshazard"` for marginal excess hazards (can only be specified if
  a relative survival model has been fitted in flexsurv).

  `"quantile"` for quantiles of the marginal all-cause survival
  distribution. The `quantiles` option also needs to be provided.

- t:

  Times to calculate marginal values at.

- ci:

  Should confidence intervals be calculated? Defaults to FALSE

- se:

  Should standard errors be calculated? Defaults to FALSE

- boot:

  Should bootstrapping be used to calculate standard error and
  confidence intervals? Defaults to FALSE, in which case the delta
  method is used

- B:

  Number of bootstrap simulations from the normal asymptotic
  distribution of the estimates used to calculate confidence intervals
  or standard errors. Decrease for greater speed at the expense of
  accuracy. Only specify if `boot = TRUE`

- cl:

  Width of symmetric confidence intervals, relative to 1.

- trans:

  Transformation to apply when calculating standard errors via the delta
  method to obtain confidence intervals. The default transformation is
  "log". Other possible names are "none", "loglog", "logit".

- contrast:

  Contrasts between standardized measures defined by `at` scenarios.
  Options are `"difference"` and `"ratio"`. There will be n-1 new
  columns created where n is the number of `at` scenarios. Default is
  NULL (i.e. no contrasts are calculated).

- trans.contrast:

  Transformation to apply when calculating standard errors for contrasts
  via the delta method to obtain confidence intervals. The default
  transformation is "none" for differences in survival, hazard,
  quantiles, or RMST, and "log" for ratios of survival, hazard,
  quantiles or RMST.

- seed:

  The random seed to use (for bootstrapping confidence intervals)

- rmap:

  An list that maps data set names to expected ratetable names. This
  must be specified if all-cause survival and hazards are required after
  fitting a relative survival model. This can also be specified if
  expected rates are required for plotting purposes. See the details
  section below.

- ratetable:

  A table of expected event rates (see
  [`ratetable`](https://rdrr.io/pkg/survival/man/ratetable.html))

- scale.ratetable:

  Transformation from the time scale of the fitted flexsurv model to the
  time scale in `ratetable`. For example, if the analysis time of the
  fitted model is in years and the ratetable is in units/day then we
  should use `scale.ratetable = 365.25`. This is the default as often
  the ratetable will be in units/day (see example).

- n.gauss.quad:

  Number of Gaussian quadrature points used for integrating the
  all-cause survival function when calculating RMST in a relative
  survival framework (default = 100)

- quantiles:

  If `type="quantile"`, this specifies the quantiles of the survival
  time distribution to return estimates for.

- interval:

  Interval of survival times for quantile root finding. Default is
  c(1e-08, 500).

## Value

A `tibble` containing one row for each time-point. The column naming
convention is `at{i}` for the ith scenario with corresponding confidence
intervals (if specified) named `at{i}_lci` and `at{i}_uci`. Contrasts
are named `contrast{k}_{j}` for the comparison of the kth versus the jth
`at` scenario.

In addition tidy long-format data.frames are returned in the attributes
`standsurv_at` and `standsurv_contrast`. These can be passed to `ggplot`
for plotting purposes (see
[`plot.standsurv`](http://chjackson.github.io/flexsurv-dev/reference/plot.standsurv.md)).

## Details

The syntax of `standsurv` follows closely that of Stata's `standsurv`
command written by Paul Lambert and Michael Crowther. The function
calculates standardized (marginal) measures including standardized
survival functions, standardized restricted mean survival times,
quantiles and the hazard of standardized survival. The standardized
survival is defined as \$\$S_s(t\|X=x) = E(S(t\|X=x,Z)) = \frac{1}{N}
\sum\_{i=1}^N S(t\|X=x,Z=z_i)\$\$ The hazard of the standardized
survival is a weighted average of individual hazard functions at time t,
weighted by the survival function at this time: \$\$h_s(t\|X=x) =
\frac{\sum\_{i=1}^N S(t\|X=x,Z=z_i)h(t\|X=x,Z=z_i)}{\sum\_{i=1}^N
S(t\|X=x,Z=z_i)}\$\$ Marginal expected survival and hazards can be
calculated by providing a population-based lifetable of class ratetable
in `ratetable` and a mapping between stratification factors in the
lifetable and the user dataset using `rmap`. If these stratification
factors are not in the fitted survival model then the user must specify
them in `newdata` along with the covariates of the model. The marginal
expected survival is calculated using the "Ederer" method that assumes
no censoring as this is most relevant approach for forecasting (see
[`survexp`](https://rdrr.io/pkg/survival/man/survexp.html)). A worked
example is given below.

Marginal all-cause survival and hazards can be calculated after fitting
a relative survival model, which utilise the expected survival from a
population ratetable. See Rutherford et al. (Chapter 6) for further
details.

## References

Paul Lambert, 2021. "STANDSURV: Stata module to compute standardized
(marginal) survival and related functions," Statistical Software
Components S458991, Boston College Department of Economics.
https://ideas.repec.org/c/boc/bocode/s458991.html

Rutherford, MJ, Lambert PC, Sweeting MJ, Pennington B, Crowther MJ,
Abrams KR, Latimer NR. 2020. "NICE DSU Technical Support Document 21:
Flexible Methods for Survival Analysis"
https://nicedsu.sites.sheffield.ac.uk/tsds/flexible-methods-for-survival-analysis-tsd

## Author

Michael Sweeting \<mikesweeting79@gmail.com\>

## Examples

``` r
## mean age is higher in those with smaller observed survival times 
newbc <- bc
set.seed(1)
newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE),
 sd = 5)

## Fit a Weibull flexsurv model with group and age as covariates
weib_age <- flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc, 
                       dist="weibull")
                       
## Calculate standardized survival and the difference in standardized survival
## for the three levels of group across a grid of survival times                        
standsurv_weib_age <- standsurv(weib_age, 
                                           at = list(list(group="Good"), 
                                                     list(group="Medium"), 
                                                     list(group="Poor")), 
                                           t=seq(0,7, length.out=100),
                                           contrast = "difference", ci=FALSE)
standsurv_weib_age
#> # A tibble: 100 × 6
#>      time   at1   at2   at3 contrast2_1 contrast3_1
#>     <dbl> <dbl> <dbl> <dbl>       <dbl>       <dbl>
#>  1 0      1     1     1         0           0      
#>  2 0.0707 0.999 0.998 0.996    -0.00111    -0.00353
#>  3 0.141  0.998 0.995 0.988    -0.00293    -0.00931
#>  4 0.212  0.996 0.991 0.980    -0.00517    -0.0164 
#>  5 0.283  0.994 0.986 0.970    -0.00772    -0.0244 
#>  6 0.354  0.992 0.981 0.959    -0.0105     -0.0331 
#>  7 0.424  0.989 0.976 0.947    -0.0135     -0.0425 
#>  8 0.495  0.987 0.970 0.935    -0.0167     -0.0523 
#>  9 0.566  0.984 0.964 0.922    -0.0201     -0.0626 
#> 10 0.636  0.981 0.958 0.908    -0.0236     -0.0732 
#> # ℹ 90 more rows

## Calculate hazard of standardized survival and the marginal hazard ratio
## for the three levels of group across a grid of survival times
## 10 bootstraps for confidence intervals (this should be larger)
if (FALSE) { # \dontrun{          
haz_standsurv_weib_age <- standsurv(weib_age, 
                                           at = list(list(group="Good"), 
                                                     list(group="Medium"), 
                                                     list(group="Poor")), 
                                           t=seq(0,7, length.out=100),
                                           type="hazard",
                                           contrast = "ratio", boot = TRUE,
                                           B=10, ci=TRUE)
haz_standsurv_weib_age                                            
plot(haz_standsurv_weib_age, ci=TRUE)
## Hazard ratio plot shows a decreasing marginal HR 
## Whereas the conditional HR is constant (model is a PH model)
plot(haz_standsurv_weib_age, contrast=TRUE, ci=TRUE)

## Calculate standardized survival from a Weibull model together with expected
## survival matching to US lifetables

# age at diagnosis in days. This is required to match to US ratetable, whose
# timescale is measured in days
newbc$agedays <- floor(newbc$age * 365.25)  
## Create some random diagnosis dates centred on 01/01/2010 with SD=1 year
## These will be used to match to expected rates in the lifetable
newbc$diag <- as.Date(floor(rnorm(dim(newbc)[1], 
                     mean = as.Date("01/01/2010", "%d/%m/%Y"), sd=365)), 
                     origin="1970-01-01")
## Create sex (assume all are female)
newbc$sex <- factor("female")
standsurv_weib_expected <- standsurv(weib_age, 
                                           at = list(list(group="Good"), 
                                                     list(group="Medium"), 
                                                     list(group="Poor")), 
                                           t=seq(0,7, length.out=100),
                                           rmap=list(sex = sex,
                                                     year = diag,
                                                     age = agedays),
                                           ratetable = survival::survexp.us,
                                           scale.ratetable = 365.25,
                                           newdata = newbc)
## Plot marginal survival with expected survival superimposed                                            
plot(standsurv_weib_expected, expected=TRUE)
} # }
```
