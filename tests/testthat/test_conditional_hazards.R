test_that("conditioning preserves hazards after the entry time", {
  d <- data.frame(time=rep(1:3, 2), x=rep(0:1, each=3))
  fit <- flexsurvreg(Surv(time) ~ x, data=d, dist="exp",
                    inits=c(.2, log(2)), fixedpars=TRUE)
  h <- summary(fit, newdata=list(x=0), type="hazard",
               t=c(4, 5, 10), start=5, ci=FALSE)[[1]]$est
  expect_equal(h, c(0, .2, .2))
  hr <- hr_flexsurvreg(fit, newdata=data.frame(x=0:1),
                      t=10, start=5, ci=FALSE)
  expect_equal(hr$est, 2)
})
