em_supported <- TRUE

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

if (em_supported){

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
  expect_equivalent(x$loglik, xd$loglik, tol=1e-06)
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
 expect_equivalent(x$loglik, xd$loglik, tol=1e-06)
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
   expect_equivalent(x$loglik, xd$loglik, tol=1e-06)
 })
 
}


test_that("Variances in Aalen-Johansen",{
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                 dists=c("gamma","gamma"),  method="em")
  aj <- ajfit_flexsurvmix(x, B=3)
  expect_true(inherits(aj, "data.frame")  & is.numeric(aj$lower))
})


