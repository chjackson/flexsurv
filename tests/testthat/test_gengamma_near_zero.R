test_that("generalized gamma is smooth near the lognormal limit", {
  for (Q in c(-1e-8, 0, 1e-8))
    expect_equal(dgengamma(c(.5, 1, 2), 0, 1, Q, log=TRUE),
                 dlnorm(c(.5, 1, 2), log=TRUE), tolerance=1e-7)
  ld <- function(Q) dgengamma(exp(2), 0, 1, Q, log=TRUE)
  expect_equal(numDeriv::grad(ld, 0), -4/3, tolerance=1e-6)
  expect_equal(numDeriv::hessian(ld, 0)[1,1], -1.5, tolerance=1e-5)
})

test_that("near-lognormal fitting retains likelihood and information", {
  set.seed(31)
  t <- rlnorm(1000)
  fit <- flexsurvreg(Surv(t) ~ 1, dist="gengamma")
  expect_equal(fit$loglik, sum(fit$logliki), tolerance=1e-8)
  expect_equal(fit$loglik, -1462.88691736, tolerance=1e-6)
  expect_gt(fit$res["Q", "se"], .05)
})
