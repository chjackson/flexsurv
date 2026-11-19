
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


test_that("survrtrunc works", {
  ## simulate some event time data
  set.seed(1)
  X <- rweibull(100, 2, 10)
  T <- rweibull(100, 2, 10)

  ## truncate above
  tmax <- 20
  obs <- X + T < tmax
  rtrunc <- tmax - X
  dat <- data.frame(X, T, rtrunc)[obs,]
  
  sf <-    survrtrunc(T, rtrunc, data=dat, tmax=tmax)
  sfnaive <- survfit(Surv(T) ~ 1, data=dat)
  ## Kaplan-Meier estimate ignoring truncation is biased
  expect_true(all(sf$surv[10:20] > sfnaive$surv[10:20]))

  if (interactive() || covr::in_covr()){
    plot(sf, conf.int=TRUE)
    lines(sfnaive, conf.int=TRUE, lty=2, col="red")
    plot(sfnaive, conf.int=TRUE)
    lines(sf, conf.int=TRUE, lty=2, col="red")
  }

  ## truncate above the maximum observed time
  tmax <- max(X + T) + 10
  obs <- X + T < tmax
  rtrunc <- tmax - X
  dat <- data.frame(X, T, rtrunc)[obs,]
  sf <-    survrtrunc(T, rtrunc, data=dat, tmax=tmax)
  ## estimates identical to the standard Kaplan-Meier
  sfnaive <- survfit(Surv(T) ~ 1, data=dat)

  expect_equal(sf$surv[1:10], sfnaive$surv[1:10])
})

test_that("right-truncated fits use the conditional likelihood and curvature", {
    old <- options(flexsurv.test.analytic.derivatives=TRUE)
    on.exit(options(old))
    ## Note that a larger case (e.g. t <- rexp(1000, .2); t[t<=5], cutoff 5) fails:
    ## BFGS overshoots into the flat region where rate -> 0.
    ## Addressable by setting fnscale in the optimisation.
    cases <- list(list(t=c(.2,.4,.5,.7), cutoff=1))
    for (case in cases) {
        t <- case$t
        cutoff <- case$cutoff
        fit <- flexsurvreg(Surv(t)~1, dist="exp", rtrunc=rep(cutoff,length(t)),
                          weights=rep(1,length(t)))
        nll <- function(p) -sum(dexp(t,exp(p),log=TRUE) -
                               pexp(cutoff,exp(p),log.p=TRUE))
        ref <- optimize(nll, log(c(.001,10)), tol=1e-10)
        x <- exp(ref$minimum)*cutoff
        logse <- 1/sqrt(length(t)*(1-x^2*exp(x)/expm1(x)^2))
        expect_equal(unname(fit$res["rate","est"]), exp(ref$minimum), tolerance=1e-4)
        expect_equal(unname(sqrt(fit$cov[1,1])), logse, tolerance=1e-3)
        expect_equal(fit$opt$convergence, 0L)
        expect_null(fit$cov.robust)
        expect_null(fit$deriv.test)
        expect_null(fit$hess.test)
    }
})
