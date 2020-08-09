## simulate events following hospital 
n <- 1000
set.seed(1)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
events <- c("icu","death","discharge")
pbase <- c(0.2, 0.3, 0.5)
event <- numeric(n)
for (i in 1:n){
  p <- pmnlogit(qmnlogit(pbase) + 2*x[i] + 3*y[i])
  event[i] <- sample(events, size=1, prob=p, replace=TRUE)
}
t <- numeric(n)
t[event=="death"] <- rgamma(sum(event=="death"), 2.5, 1.2)
t[event=="discharge"] <- rgamma(sum(event=="discharge"), 3.5, 0.6)
t[event=="icu"] <- rgamma(sum(event=="icu"), 1, 3.2)
cens <- as.numeric(t > 3)
t[cens] <- 3
status <- 1 - cens
dat <- data.frame(t, status, x, y, event)

## model for event following hospital
fhosp <-   flexsurvmix(Surv(t, status) ~ x, pformula = ~x + y,
                       data=dat, event=event, 
                       dists=c("gamma","gamma","gamma"))


## simulate events following ICU
nicu <- sum(dat$event=="icu")
picu <- c(0.4, 0.6)
set.seed(1)
evicu <- sample(c("death","discharge"), size=nicu, prob=picu, replace=TRUE)
ti <- numeric(nicu)
ti[evicu=="death"] <- rgamma(sum(evicu=="death"), 1.5, 1)
ti[evicu=="discharge"] <- rgamma(sum(evicu=="discharge"), 0.5, 3)
censi <- as.numeric(ti > 1)
ti[censi] <- 1
statusi <- 1 - censi
dati <- data.frame(ti, statusi, evicu)

## model for event following ICU
ficu <- flexsurvmix(Surv(ti, statusi) ~ 1, data=dati, event=evicu, 
                       dists=c("gamma","gamma"))

## Construct multi-state model object
fm <- fmixmsm("hospital"=fhosp, "icu"=ficu)

test_that("prob_pathway",{
  probh <- probs_flexsurvmix(fhosp)
  probi <- probs_flexsurvmix(ficu)
  pp <- prob_pathway(fm)
  expect_equal(
    probh$val[probh$event=="icu"] * probi$val[probi$event=="death"],
    pp$val[pp$pathway=="hospital-icu-death"]
  )
  prob_pathway(fm, final=TRUE)
  nd <- data.frame(x=c(0,0.02), y=c(0,0.01))
  probh <- probs_flexsurvmix(fhosp,newdata=nd)
  probi <- probs_flexsurvmix(ficu,newdata=nd)
  pp <- prob_pathway(fm, newdata=nd)
  expect_equal(
    probh$val[probh$event=="icu" & probh$x==0.02 & probh$y==0.01] * 
      probi$val[probi$event=="death" & probi$x==0.02 & probi$y==0.01],
    pp$val[pp$pathway=="hospital-icu-death" & pp$x==0.02 & pp$y==0.01]
  )
  prob_pathway(fm, newdata=nd, final=TRUE)
  expect_true(is.numeric(prob_pathway(fm, B=3)$lower))
  expect_true(is.numeric(prob_pathway(fm, final=TRUE, B=3)$lower))
  expect_true(is.numeric(prob_pathway(fm, newdata=nd, B=3)$lower))
  expect_true(is.numeric(prob_pathway(fm, newdata=nd, final=TRUE, B=3)$lower))
})

test_that("mean_tofinal",{
  mean_tofinal(fm)
  mean_tofinal(fm, final=TRUE)
  mean_tofinal(fm, newdata=nd)
  mean_tofinal(fm, newdata=nd, final=TRUE)
  expect_true(is.numeric(mean_tofinal(fm, B=3)$lower))
  expect_true(is.numeric(mean_tofinal(fm, final=TRUE, B=3)$lower))
  expect_true(is.numeric(mean_tofinal(fm, newdata=nd, B=3)$lower))
  expect_true(is.numeric(mean_tofinal(fm, newdata=nd, final=TRUE, B=3)$lower))
})