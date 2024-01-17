## Thanks to Andrea Discacciati @https://github.com/anddis
## Fixed https://github.com/chjackson/flexsurv/issues/178

test_that("Non-default factor contrasts", { 
  veteran$celltype <- factor(veteran$celltype)
  fit.noc <- flexsurvreg(Surv(time, status) ~ celltype,
                         data = veteran, 
                         dist = "exp")

  veteran$celltype2 <- veteran$celltype
  contrasts(veteran$celltype2) <- stats::contr.helmert(4)
  fit.c <- flexsurvreg(Surv(time, status) ~ celltype2,
                       data = veteran, 
                       dist = "exp")

  fit.poi <- glm(status ~ celltype2 + offset(log(time)),
                 data = veteran,
                 family = "poisson")

  expect_equivalent(coef(fit.c), coef(fit.poi))
  expect_true(coef(fit.c)[1] != coef(fit.noc)[1])

  summ.noc <- summary(fit.noc, type="survival", t=10, tidy=TRUE, ci=FALSE)
  summ.c <- summary(fit.c, type="survival", t=10, tidy=TRUE, ci=FALSE)
  expect_equivalent(summ.noc$estimates, summ.c$estimates)
})
