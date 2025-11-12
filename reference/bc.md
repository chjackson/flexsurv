# Breast cancer survival data

Survival times of 686 patients with primary node positive breast cancer.

## Usage

``` r
bc
```

## Format

A data frame with 686 rows.

|           |           |                                                                    |
|-----------|-----------|--------------------------------------------------------------------|
| `censrec` | (numeric) | 1=dead, 0=censored                                                 |
| `rectime` | (numeric) | Time of death or censoring in days                                 |
| `group`   | (numeric) | Prognostic group: `"Good"`,`"Medium"` or `"Poor"`,                 |
|           |           | from a regression model developed by Sauerbrei and Royston (1999). |

## Source

German Breast Cancer Study Group, 1984-1989. Used as a reference dataset
for the spline-based survival model of Royston and Parmar (2002),
implemented here in
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).
Originally provided with the `stpm` (Royston 2001, 2004) and `stpm2`
(Lambert 2009, 2010) Stata modules.

## References

Royston, P. and Parmar, M. (2002). Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of
treatment effects. Statistics in Medicine 21(1):2175-2197.

Sauerbrei, W. and Royston, P. (1999). Building multivariable prognostic
and diagnostic models: transformation of the predictors using fractional
polynomials. Journal of the Royal Statistical Society, Series A
162:71-94.

## See also

[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
