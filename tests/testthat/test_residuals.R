test_that("residuals",{
  fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
  expect_true(is.numeric(residuals(fitg, type="response")))
  expect_true(is.numeric(residuals(fitg, type="coxsnell")))
  cs <- coxsnell_flexsurvreg(fitg)
  surv <- survfit(Surv(cs$est, ovarian$fustat) ~ 1)
  if (interactive()){
    plot(surv, fun="cumhaz", xlim=c(0,1), ylim=c(0,1))
    abline(0, 1, col="red")
  }
  expect_lt(max(surv$cumhaz[-1] / surv$time[-1]), 1.5)
  expect_gt(min(surv$cumhaz[-1] / surv$time[-1]), 0.5)
})
