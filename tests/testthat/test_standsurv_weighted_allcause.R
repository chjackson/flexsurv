test_that("weighted all-cause hazards add excess and background hazards", {
  d <- data.frame(time=1:3, bh=.1)
  fit <- flexsurvreg(Surv(time) ~ 1, data=d, dist="exp", bhazard=bh,
                    inits=.2, fixedpars=TRUE)
  nd <- data.frame(age=rep(65*365.25, 2), sex=factor(rep("female", 2)),
                   year=as.Date(rep("2010-01-01", 2)))
  nd[["(weights)"]] <- c(1, 2)
  rt <- survival::survexp.us
  rt[] <- .1/365.25
  ss <- suppressMessages(standsurv(fit, newdata=nd, t=1, type="hazard",
      rmap=list(age=age, sex=sex, year=year), ratetable=rt, scale.ratetable=365.25))
  expect_equal(ss$at1, .3, tolerance=1e-6)
})
