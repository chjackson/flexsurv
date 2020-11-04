
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
t[cens] <- 4 
dat <- data.frame(t, status, event)
dat$evname <- c("cure", "death")[dat$event]
dat$evnamef <- factor(dat$evname)

test_that("flexsurvmix basic",{
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"))
  expect_equivalent(fs$loglik, -1343.37172263181)
  
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), fixedpars=TRUE)
  expect_equivalent(fs$loglik, -1550.65934372248)
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
t[cens] <- 2
status <- 1 - cens
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x)
dat$evname <- c("cure", "death")[dat$event]


test_that("flexsurvmix with covariates on time to event distributions",{
  
  ## covariates on all location pars
  fs <- flexsurvmix(Surv(t, status) ~ x, data=dat, event=event, 
                    dists=c("gamma","weibull"), fixedpars=TRUE)
  
  ## covariates on all location pars, and some shape pars
  ## this is the right model here 
  fs <- flexsurvmix(Surv(t, status) ~ x, data=dat, event=event, dists=c("gamma","gamma"), 
                    anc=list(list(shape=~x), list(shape=~1)), fixedpars=TRUE)
  
  ## covariates supplied manually for location pars
  fs <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), 
                    anc=list(list(rate=~x, shape=~x), list(rate=~x)), fixedpars=TRUE)
  
  ## covariates on location for one component but not another. two ways to do
  ini <- list(c(1,3,0.6,0.2), c(1,0.3))
  fs1 <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=event, dists=c("gamma","gamma"), 
                     anc=list(list(rate=~x, shape=~x), list(rate=~1)), inits=ini, fixedpars=TRUE)
  
  fs2 <- flexsurvmix(list(`1`=Surv(t, status) ~ x, 
                          `2`=Surv(t, status) ~ 1),
                     data=dat, event=event, dists=c("gamma","gamma"), 
                     anc=list(list(shape=~x), list(rate=~1)), inits=ini, fixedpars=TRUE)
  expect_equal(fs1$loglik, fs2$loglik)
  
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
t[cens] <- 3
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x)
dat$evname <- c("cure", "death")[dat$event]

test_that("flexsurvmix with covariates on mixing probabilities",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), fixedpars=TRUE)
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), pformula = ~ x + y, fixedpars=TRUE)
})


## Partially-observed outcomes
n <- 1000
set.seed(1)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
events <- c("icu","death","discharge")
p <- c(0.2, 0.3, 0.5)
event <- sample(events, size=n, prob=p, replace=TRUE)
t <- numeric(n)
t[event=="death"] <- rgamma(sum(event=="death"), 2.5, 1.2)
t[event=="discharge"] <- rgamma(sum(event=="discharge"), 3.5, 0.6)
t[event=="icu"] <- rgamma(sum(event=="icu"), 1, 3.2)
cens <- as.numeric(t > 3)
t[cens] <- 3
status <- 1 - cens
##  Out of those still alive at cens time, or those discharged, label some as as partially observed, 
## so we don't know if they've been discharged yet or they are still in hosp at the cens time
pobs_eligible <- (event == "discharge") | (t > 3)
pobs <- rbinom(n=sum(pobs_eligible), size=1, prob=0.2) == 1
## Construct data  
# Partially-obs people: know alive at cens time, else don't know anything
# Time of death and ICU are right cens 
# Time of discharge interval cens on 0 to Inf ie no info 
# Note we don't need a "partial_events" indicator in the event data 
# partially obs people still at risk of all events in the future 
# we just know they are still alive at cens time, so could die in future or go to ICU 
t1disc <- t2disc <- t
t2disc[cens] <- Inf
t1disc[pobs] <- 0
dat <- data.frame(t, status, t1disc, t2disc, event) 
