test_that("bootstrap contrasts preserve shared parameter draws", {
  d <- data.frame(time=c((1:50)/10, (1:50)/20), x=rep(0:1, each=50))
  fit <- flexsurvreg(Surv(time) ~ x, data=d, dist="exp")
  B <- 200
  ss <- suppressMessages(standsurv(fit, t=2,
      at=list(list(x=0), list(x=.25)), contrast="difference",
      boot=TRUE, B=B, seed=11, ci=TRUE, se=TRUE))
  set.seed(11)
  sim <- normboot.flexsurvreg(fit, B=B, newdata=data.frame(x=c(0, .25)))
  difference <- exp(-2*sim[[2]][,"rate"]) - exp(-2*sim[[1]][,"rate"])
  expect_equal(ss$contrast2_1_se, sd(difference))
  expect_equal(c(ss$contrast2_1_lci, ss$contrast2_1_uci),
               unname(quantile(difference, c(.025, .975))))
})
