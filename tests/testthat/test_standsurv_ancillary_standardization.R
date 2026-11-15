test_that("standardization includes ancillary-only covariates", {
  d <- data.frame(time=1:4, age=c(0, 0, 1, 1))
  fit <- flexsurvreg(Surv(time) ~ 1, anc=list(shape=~age),
                    data=d, dist="weibull", inits=c(1, 1, log(4)),
                    fixedpars=TRUE)
  nd <- data.frame(age=0:1)
  expected <- mean(pweibull(.5, shape=c(1, 4), scale=1, lower.tail=FALSE))
  ss <- standsurv(fit, newdata=nd, t=.5)
  reversed <- standsurv(fit, newdata=nd[2:1,,drop=FALSE], t=.5)
  expect_equal(ss$at1, expected)
  expect_equal(reversed$at1, expected)
})
