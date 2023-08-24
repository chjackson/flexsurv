if (!identical(Sys.getenv("NOT_CRAN"), "true")) return()

n <- 200
p <- 0.5 
set.seed(1)
death <- rbinom(n, 1, p)
t <- numeric(n)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
t[death==0] <- rgamma(sum(death==0), 1*exp(0.3*x[death==0]), 3.2*exp(0.5*x[death==0]))
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
cens <- as.numeric(t > 2)
t[cens] <- 2
status <- 1 - cens
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x)
dat$evname <- c("cure", "death")[dat$event]

test_that("EM basic",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"),  method="em")
  xd <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                    dists=c("gamma","gamma"),  method="direct")
  expect_equivalent(x$loglik, xd$loglik, tol=1e-06)
})

test_that("EM pformula",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), pformula = ~x+y,  method="em")
  xd <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                    dists=c("gamma","gamma"), pformula = ~x+y,  method="direct")
  expect_equivalent(x$loglik, xd$loglik, tol=1e-04)
})

test_that("EM with different covs on different components", {
  x <- flexsurvmix(list(cure=Surv(t, status) ~ 1, 
                        death=Surv(t, status) ~ x),
                   data=dat, event=evname, 
                   dists=c("gamma","gamma"),  method="em")
  xd <- flexsurvmix(list(cure=Surv(t, status) ~ 1, 
                         death=Surv(t, status) ~ x),
                    data=dat, event=evname, 
                    dists=c("gamma","gamma"),  method="direct")
  expect_equivalent(x$loglik, xd$loglik, tol=1e-04)
})

test_that("EM with anc and different covs on different components", {
  x <- flexsurvmix(list(cure=Surv(t, status) ~ 1, 
                        death=Surv(t, status) ~ x),
                   data=dat, event=evname, 
                   anc=list(cure=list(shape=~y), death=list(shape=~1)),
                   dists=c("gamma","gamma"),  method="em")
  xd <- flexsurvmix(list(cure=Surv(t, status) ~ 1, 
                         death=Surv(t, status) ~ x),
                    data=dat, event=evname, 
                    anc=list(cure=list(shape=~y), death=list(shape=~1)),
                    dists=c("gamma","gamma"),  method="direct")
  expect_equivalent(x$loglik, xd$loglik, tol=1e-04)
})

test_that("Variances in Aalen-Johansen",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"),  method="em")
  aj <- ajfit_flexsurvmix(x, B=3)
  expect_true(inherits(aj, "data.frame")  & is.numeric(aj$lower))
})

n <- 100
set.seed(1)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
p <- plogis(qlogis(0.5) + 2*x - 3*y)
death <- rbinom(n, 1, p)
t <- numeric(n)
t[death==0] <- rgamma(sum(death==0), 1, 3.2)
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
cens <- as.numeric(t > 3)
status <- 1 - cens
t[cens] <- 3
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x)
dat$evname <- c("cure", "death")[dat$event]

test_that("var.method",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), pformula = ~ x + y, method="em", fixedpars=c(3,4),
                   em.control=list(var.method="louis"))
  x
  
  x2 <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), pformula = ~ x + y, method="em", fixedpars=c(3,4),
                   em.control=list(var.method="direct"))
  
  x2
  expect_equivalent(x$loglik, x2$loglik, tol=1e-04)
  expect_true(x$cov[1,1] != x2$cov[1,1])
})

test_that("Models with interactions",{
  dat$x <- factor(rbinom(nrow(dat), size=1, prob=0.5))
  dat$y <- factor(rbinom(nrow(dat), size=1, prob=0.4))
  fsi <- flexsurvmix(formula = Surv(t, status) ~ x*y, data = dat, event = event, 
                     dists = c("gamma", "weibull"), fixedpars = FALSE, method="em", 
                     em.control=list(var.method="louis", reltol=1e-04))
  aj <- ajfit_flexsurvmix(fsi)
  expect_equal(aj$val[1], 1)  
})
