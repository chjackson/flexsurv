test_that("Gompertz analytical derivatives and Hessian match numeric values at and near shape = 0", {
  t <- c(0.5, 2, 10)
  rate <- 0.1
  for (sh in c(-1e-6, -1e-8, 0, 1e-8, 1e-6, 0.05, 0.5)) {
    ld <- function(p) sum(dgompertz(t, shape = p[1], rate = exp(p[2]), log = TRUE))
    ls <- function(p) sum(pgompertz(t, shape = p[1], rate = exp(p[2]), lower.tail = FALSE, log.p = TRUE))
    p0 <- c(sh, log(rate))

    an_d1 <- colSums(flexsurv:::DLdgompertz(t, rep(sh, length(t)), rep(rate, length(t))))
    an_s1 <- colSums(flexsurv:::DLSgompertz(t, rep(sh, length(t)), rep(rate, length(t))))
    expect_equal(unname(an_d1), numDeriv::grad(ld, p0), tolerance = 1e-6)
    expect_equal(unname(an_s1), numDeriv::grad(ls, p0), tolerance = 1e-6)

    an_d2 <- apply(flexsurv:::D2Ldgompertz(t, rep(sh, length(t)), rep(rate, length(t))), 2:3, sum)
    an_s2 <- apply(flexsurv:::D2LSgompertz(t, rep(sh, length(t)), rep(rate, length(t))), 2:3, sum)
    expect_equal(unname(an_d2), numDeriv::hessian(ld, p0), tolerance = 1e-5)
    expect_equal(unname(an_s2), numDeriv::hessian(ls, p0), tolerance = 1e-5)
  }
})

test_that("Gompertz fit starting at shape = 0 converges to the same MLE and SEs", {
  fit_zero <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "gompertz",
                          inits = c(0, 0.001))
  fit_pos <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "gompertz",
                         inits = c(0.001, 0.001))
  expect_equal(fit_zero$res[, "est"], fit_pos$res[, "est"], tolerance = 1e-4)
  expect_equal(fit_zero$res[, "se"], fit_pos$res[, "se"], tolerance = 1e-3)
  expect_equal(fit_zero$loglik, fit_pos$loglik, tolerance = 1e-6)
})
