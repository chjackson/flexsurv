# helper function to safely convert a Hessian matrix to covariance matrix

helper function to safely convert a Hessian matrix to covariance matrix

## Usage

``` r
.hess_to_cov(hessian, tol.solve = 1e-09, tol.evalues = 1e-05, ...)
```

## Arguments

- hessian:

  hessian matrix to convert to covariance matrix (must be evaluated at
  MLE)

- tol.solve:

  tolerance used for solve()

- tol.evalues:

  accepted tolerance for negative eigenvalues of the covariance matrix

- ...:

  arguments passed to Matrix::nearPD
