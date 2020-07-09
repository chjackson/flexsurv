
set.seed(1)
## simulate time to initial event
X <- rexp(1000, 0.2)
## simulate time between initial and final event
tdelay <- rgamma(1000, 2, 10)

tmax <- 40
obs <- X + tdelay < tmax 
rtrunc <- tmax - X
dat <- data.frame(X, tdelay, rtrunc)[obs,]

fs <- flexsurvrtrunc(t=tdelay, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
                     dist="gamma", theta=0.2)

test_that("flexsurvrtrunc works", {
  expect_equal(fs$loglik, -7846.25011895056)
})

test_that("summary.flexsurvtrunc works", { 
  summmean <- summary(fs, type="mean")
  expect_equal(summmean$est, 0.2095068)
  summmed <- summary(fs, type="median")
  expect_equal(summmed$est, 0.1746857)
  summary(fs, type="rmst", t=0.3)
  summq <- summary(fs, type="quantile")
  summq2 <- summary(fs, type="quantile", quantiles=c(0.025, 0.5, 0.975))
  expect_equal(summmed$est, summq$est)
  expect_equal(summmed$est, summq2$est[summq2$quantile==0.5])
  summ <-  summary(fs, type="survival", t=c(0.01, 0.02, 0.03))
  expect_equal(summ$est[summ$time==0.02], 
               pgamma(0.02, fs$res["shape","est"], fs$res["rate","est"], lower.tail=FALSE))
  summary(fs, type="cumhaz", t=seq(0.01, 0.05, by=0.01))
  summary(fs, type="hazard", t=seq(0.01, 0.05, by=0.01))
  
  fntest <- function(shape, rate){2 * mean_gamma(shape,rate)}
  summfn <- summary(fs, fn=fntest, t=1)
  expect_equal(summfn$est, 2*summmean$est)
  
  set.seed(1)
  summse <- summary(fs, type="median", se=TRUE)
  expect_equal(summse$se, 0.004329825)
})

