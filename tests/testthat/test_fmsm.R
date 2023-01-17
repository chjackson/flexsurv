#### Functions that summarise outputs from fmsm multi-state model objects

## Three state competing risks model for testing functions that only work on
## Markov models

tmat <- rbind(c(NA,1,2),c(NA,NA,NA),c(NA,NA,NA))
set.seed(1)
bosms3$x <- rnorm(nrow(bosms3))

bweic <- bweim <- vector(2, mode="list")
for (i in 1:2) {
  bweic[[i]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==i),
                            data = bosms3, dist = "weibull")
  bweim[[i]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==i),
                            data = bosms3, dist = "weibull")
}

weic <- fmsm("Well-BOS"=bweic[[1]], "Well-Death"=bweic[[2]], trans=tmat)
weim <- fmsm("Well-BOS"=bweim[[1]], "Well-Death"=bweim[[2]], trans=tmat)
nd <- data.frame(x=c(0,0.01,-10,10))

test_that("pfinal_fmsm", { 
  expect_equal(pfinal_fmsm(weim, fromstate="State 1")$val, c(0.717737627151196, 0.282262372848804))
  expect_error(pfinal_fmsm(weim, fromstate="State 2"), "No destination states")
  expect_equal(pfinal_fmsm(weic, newdata=nd, fromstate="State 1")$val[1:2], c(0.715828147070156, 0.715333178039627))
  expect_true(is.numeric(pfinal_fmsm(weim, fromstate="State 1", B=3)$lower))
  expect_equal(pfinal_fmsm(weim, fromstate="State 1", maxt=100000)$val,
               pfinal_fmsm(weim, fromstate="State 1", maxt=10000000)$val, tolerance=1e-06)
  expect_error(pfinal_fmsm(weim, fromstate="1"), "not found")
})

test_that("simfinal_fmsm",{
  set.seed(1)
  sm <- simfinal_fmsm(weim)
  expect_equal(sm$val[sm$quantity=="prob"], c(0.71758, 0.28242))
  expect_equal(sm$val[sm$quantity=="50%"], c(2.25618202292384, 4.67377970752349))
  sm2 <- simfinal_fmsm(weim, probs=c(0.25, 0.75))
  sm3 <- simfinal_fmsm(weim, probs=c(0.25, 0.75), t=10000)
  expect_equal(sm2$val[sm$quantity=="prob"], sm3$val[sm$quantity=="prob"], tolerance=0.1)
  sm2 <- simfinal_fmsm(weim, probs=c(0.25, 0.75), M=1000, B=10)
  
  expect_error(simfinal_fmsm(weic), "`newdata` should be supplied")
  nd <- data.frame(x=c(0, 0.01))
  simfinal_fmsm(weic, newdata=nd)
  simfinal_fmsm(weic, newdata=nd, M=1000, B=10)
})
