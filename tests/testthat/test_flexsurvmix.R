
n <- 1000
p <- 0.5 
set.seed(1)
death <- rbinom(n, 1, p)
t <- numeric(n)
t[death==0] <- rgamma(sum(death==0), 1, 3.2)
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
cens <- as.numeric(t > 4)
status <- 1 - cens
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event)
dat$evname <- c("cure", "death")[dat$event]
dat$evnamef <- factor(dat$evname)

test_that("flexsurvmix basic",{
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"))
  expect_equivalent(fs$loglik, -1331.6712793566)
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), fixedpars=TRUE)
  expect_equivalent(fs$loglik, -1548.52559021999)
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), fixedpars=1:2)
  
  expect_silent({
  ## event as character
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, dists=c("gamma","gamma"), fixedpars=TRUE)
  ## event as factor
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evnamef, dists=c("gamma","gamma"), fixedpars=TRUE)
  ## user supplied inits 
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, dists=c("gamma","gamma"), 
                    inits = list(cure=c(2, 1), death=c(1, 3)), fixedpars=TRUE)
  
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","weibull"), fixedpars=TRUE)
  })
})

n <- 1000
p <- 0.5 
set.seed(1)
death <- rbinom(n, 1, p)
t <- numeric(n)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
t[death==0] <- rgamma(sum(death==0), 1*exp(0.3*x[death==0]), 3.2*exp(0.5*x[death==0]))
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
cens <- as.numeric(t > 2)
status <- 1 - cens
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x)
dat$evname <- c("cure", "death")[dat$event]

test_that("flexsurvmix with covariates on time to event distributions",{

  ## covariates on all location pars
  fs <- flexsurvmix(Surv(t, status) ~ x, data=dat, event=event, dists=c("gamma","weibull"), fixedpars=TRUE)
  
  ## covariates on all location pars, and some shape pars
  ## this is the right model here 
  fs <- flexsurvmix(Surv(t, status) ~ x, data=dat, event=event, dists=c("gamma","gamma"), 
                    anc=list(list(shape=~x), list(shape=~1)), fixedpars=TRUE)
  
  ## covariates supplied manually for location pars
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), 
                    anc=list(list(rate=~x, shape=~x), list(rate=~x)), fixedpars=TRUE)
  
  ## covariates on one component but not another 
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), 
                    anc=list(list(rate=~x, shape=~x), list(rate=~1)), fixedpars=TRUE)
  
  expect_error(
    {fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), 
                       anc=list(list(rate=~x, shape=~x)), fixedpars=TRUE)},
    "`anc` of length 1, should equal the number of")
  
  ## user-supplied inits
  fs <- flexsurvmix(Surv(t, status) ~ x, data=dat, event=evname, 
                    dists=c("gamma","gamma"), 
                    anc=list(list(shape=~x), list(shape=~1)),
                    inits=list(cure=c(2, 1, 0, 0), death=c(1, 3, 0)), fixedpars=TRUE)
  
})

n <- 10000
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
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x)
dat$evname <- c("cure", "death")[dat$event]

##  If censoring time is low enough, then get TTE estimates for death wrong. overestimated. 
## if leave out pformula, then tte ests are OK. probs misspecified 
## likely component membership determined by covariates as well as time 
## looks like identifiability thing

test_that("flexsurvmix with covariates on mixing probabilities",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), fixedpars=TRUE)
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), pformula = ~ x + y, fixedpars=TRUE)
})

