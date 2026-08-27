# Simulate from the asymptotic normal distribution of parameter estimates.

Produce a matrix of alternative parameter estimates under sampling
uncertainty, at covariate values supplied by the user. Used by
[`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md)
for obtaining confidence intervals around functions of parameters.

## Usage

``` r
normboot.flexsurvreg(
  x,
  B,
  newdata = NULL,
  X = NULL,
  transform = FALSE,
  raw = FALSE,
  tidy = FALSE,
  rawsim = NULL
)
```

## Arguments

- x:

  A fitted model from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  (or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)).

- B:

  Number of samples.

- newdata:

  Data frame or list containing the covariate values to evaluate the
  parameters at. If there are covariates in the model, at least one of
  `newdata` or `X` must be supplied, unless `raw=TRUE`.

- X:

  Alternative (less convenient) format for covariate values: a matrix
  with one row, with one column for each covariate or factor contrast.
  Formed from all the "model matrices", one for each named parameter of
  the distribution, with intercepts excluded, `cbind`ed together.

- transform:

  `TRUE` if the results should be transformed to the real-line scale,
  typically by log if the parameter is defined as positive. The default
  `FALSE` returns parameters on the natural scale.

- raw:

  Return samples of the baseline parameters and the covariate effects,
  rather than the default of adjusting the baseline parameters for
  covariates.

- tidy:

  If `FALSE` (the default) then a list is returned. If `TRUE` a data
  frame is returned, consisting of the list elements `rbind`ed together,
  with integer variables labelling the covariate number and simulation
  replicate number.

- rawsim:

  allows input of raw samples from a previous run of
  `normboot.flexsurvreg`. This is useful if running
  `normboot.flexsurvreg` multiple time on the same dataset but with
  counterfactual contrasts, e.g. treat =0 vs. treat =1. Used in
  `standsurv.flexsurvreg`.

## Value

If `newdata` includes only one covariate combination, a matrix will be
returned with `B` rows, and one column for each named parameter of the
survival distribution.

If more than one covariate combination is requested (e.g. `newdata` is a
data frame with more than one row), then a list of matrices will be
returned, one for each covariate combination.

## References

Mandel, M. (2013). "Simulation based confidence intervals for functions
with complicated derivatives." The American Statistician (in press).

## See also

[`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md)

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

## Examples

``` r

    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
#>               rate
#>  [1,] 0.0001333406
#>  [2,] 0.0002711950
#>  [3,] 0.0003700718
#>  [4,] 0.0002548242
#>  [5,] 0.0005457485
#>  [6,] 0.0002299253
#>  [7,] 0.0004888764
#>  [8,] 0.0001547695
#>  [9,] 0.0005135426
#> [10,] 0.0004136617
#> attr(,"X")
#>   age
#> 1  50
#> attr(,"X")attr(,"newdata")
#>   age
#> 1  50
#> attr(,"rawsim")
#>            rate        age
#>  [1,] -18.29776 0.18750322
#>  [2,] -14.33221 0.12239078
#>  [3,] -14.75103 0.13698433
#>  [4,] -14.90062 0.13251367
#>  [5,] -10.99314 0.06959569
#>  [6,] -15.71307 0.14670627
#>  [7,] -13.66212 0.12077431
#>  [8,] -17.76240 0.17977660
#>  [9,] -11.44435 0.07740335
#> [10,] -12.31628 0.09051640
    normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
#>               rate
#>  [1,] 0.0003826964
#>  [2,] 0.0007727804
#>  [3,] 0.0002973256
#>  [4,] 0.0004210293
#>  [5,] 0.0003696239
#>  [6,] 0.0003892618
#>  [7,] 0.0002302111
#>  [8,] 0.0004541271
#>  [9,] 0.0002247796
#> [10,] 0.0002689798
#> attr(,"X")
#>      [,1]
#> [1,]   50
#> attr(,"rawsim")
#>            rate        age
#>  [1,] -12.31256 0.08888573
#>  [2,] -10.19412 0.06057214
#>  [3,] -12.99376 0.09746149
#>  [4,] -14.20243 0.12859241
#>  [5,] -12.97912 0.10152183
#>  [6,] -13.83425 0.11965986
#>  [7,] -14.85106 0.12949092
#>  [8,] -10.52134 0.05648414
#>  [9,] -15.71354 0.14626291
#> [10,] -14.89006 0.13338376
    normboot.flexsurvreg(fite, B=10, newdata=list(age=0))  ## closer to...
#>               rate
#>  [1,] 9.179604e-06
#>  [2,] 8.001138e-06
#>  [3,] 3.313562e-06
#>  [4,] 5.777548e-08
#>  [5,] 2.844863e-06
#>  [6,] 6.954191e-07
#>  [7,] 6.404087e-07
#>  [8,] 6.187050e-06
#>  [9,] 2.166374e-06
#> [10,] 2.559143e-07
#> attr(,"X")
#>   age
#> 1   0
#> attr(,"X")attr(,"newdata")
#>   age
#> 1   0
#> attr(,"rawsim")
#>            rate        age
#>  [1,] -11.59853 0.07562969
#>  [2,] -11.73593 0.09619946
#>  [3,] -12.61749 0.09538749
#>  [4,] -16.66670 0.16569939
#>  [5,] -12.77000 0.10285181
#>  [6,] -14.17875 0.12219826
#>  [7,] -14.26116 0.12455604
#>  [8,] -11.99305 0.08907006
#>  [9,] -13.04246 0.09910768
#> [10,] -15.17842 0.13991248
    fite$res
#>               est         L95%         U95%           se
#> rate 8.883706e-07 1.440514e-08 5.478617e-05 1.868243e-06
#> age  1.185227e-01 5.210364e-02 1.849418e-01 3.388790e-02
```
