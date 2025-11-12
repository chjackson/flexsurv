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
#>  [1,] 0.0003700718
#>  [2,] 0.0002548242
#>  [3,] 0.0005457485
#>  [4,] 0.0002299253
#>  [5,] 0.0004888764
#>  [6,] 0.0001547695
#>  [7,] 0.0005135426
#>  [8,] 0.0004136617
#>  [9,] 0.0003826964
#> [10,] 0.0007727804
#> attr(,"X")
#>   age
#> 1  50
#> attr(,"X")attr(,"newdata")
#>   age
#> 1  50
#> attr(,"rawsim")
#>            rate        age
#>  [1,] -14.75103 0.13698433
#>  [2,] -14.90062 0.13251367
#>  [3,] -10.99314 0.06959569
#>  [4,] -15.71307 0.14670627
#>  [5,] -13.66212 0.12077431
#>  [6,] -17.76240 0.17977660
#>  [7,] -11.44435 0.07740335
#>  [8,] -12.31628 0.09051640
#>  [9,] -12.31256 0.08888573
#> [10,] -10.19412 0.06057214
    normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
#>               rate
#>  [1,] 0.0002973256
#>  [2,] 0.0004210293
#>  [3,] 0.0003696239
#>  [4,] 0.0003892618
#>  [5,] 0.0002302111
#>  [6,] 0.0004541271
#>  [7,] 0.0002247796
#>  [8,] 0.0002689798
#>  [9,] 0.0004028114
#> [10,] 0.0009819658
#> attr(,"X")
#>      [,1]
#> [1,]   50
#> attr(,"rawsim")
#>            rate        age
#>  [1,] -12.99376 0.09746149
#>  [2,] -14.20243 0.12859241
#>  [3,] -12.97912 0.10152183
#>  [4,] -13.83425 0.11965986
#>  [5,] -14.85106 0.12949092
#>  [6,] -10.52134 0.05648414
#>  [7,] -15.71354 0.14626291
#>  [8,] -14.89006 0.13338376
#>  [9,] -11.59853 0.07562969
#> [10,] -11.73593 0.09619946
    normboot.flexsurvreg(fite, B=10, newdata=list(age=0))  ## closer to...
#>               rate
#>  [1,] 3.313562e-06
#>  [2,] 5.777548e-08
#>  [3,] 2.844863e-06
#>  [4,] 6.954191e-07
#>  [5,] 6.404087e-07
#>  [6,] 6.187050e-06
#>  [7,] 2.166374e-06
#>  [8,] 2.559143e-07
#>  [9,] 2.512828e-06
#> [10,] 4.388557e-07
#> attr(,"X")
#>   age
#> 1   0
#> attr(,"X")attr(,"newdata")
#>   age
#> 1   0
#> attr(,"rawsim")
#>            rate        age
#>  [1,] -12.61749 0.09538749
#>  [2,] -16.66670 0.16569939
#>  [3,] -12.77000 0.10285181
#>  [4,] -14.17875 0.12219826
#>  [5,] -14.26116 0.12455604
#>  [6,] -11.99305 0.08907006
#>  [7,] -13.04246 0.09910768
#>  [8,] -15.17842 0.13991248
#>  [9,] -12.89410 0.10026367
#> [10,] -14.63910 0.12846751
    fite$res
#>               est         L95%         U95%           se    se.robust
#> rate 8.883706e-07 1.440514e-08 5.478617e-05 1.868243e-06 1.714961e-06
#> age  1.185227e-01 5.210364e-02 1.849418e-01 3.388790e-02 2.978249e-02
```
