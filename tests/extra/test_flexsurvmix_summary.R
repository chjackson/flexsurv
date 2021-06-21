

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

fs <- flexsurvmix(Surv(t, status) ~ x, data=dat, event=event, 
                  dists=c("gamma","weibull"), fixedpars=FALSE,
                  optim.control=list(reltol=1e-4))

test_that("summary functions, default newdata, one covariate",{
  expect_equal(mean_flexsurvmix(fs)$val[1], 0.337598824779061, tol=1e-06)
  qu <- quantile_flexsurvmix(fs)
  expect_equal(qu$val[qu$event==1 & qu$p==0.975], 1.262683082, tol=1e-06)
  pr <- probs_flexsurvmix(fs)
  expect_equal(pr$val[1], 0.522480219512026)
  pr <- p_flexsurvmix(fs)
  expect_equal(pr$val[pr$state=="1"], 0.494566645108882)
  pdf_flexsurvmix(fs, t=1:4)
  
  
})


test_that("summary functions, newdata, one covariate",{
  nd <- list(x=1)
  expect_equal(mean_flexsurvmix(fs,newdata=nd)$val[1], 0.280035901598242, tol=1e-06)
  qu <- quantile_flexsurvmix(fs, newdata=nd)
  expect_equal(qu$val[qu$event==1 & qu$p==0.975], 1.04738692622352, tol=1e-06)
  pr <- probs_flexsurvmix(fs, newdata=nd)
  expect_equal(pr$val[1], 0.522480219512026)
  pdf_flexsurvmix(fs, newdata=list(x=1), t=1:4)
  
  nd <- list(x=c(1,2))
  pr <- p_flexsurvmix(fs,newdata=nd)
  expect_equal(pr$val[pr$state=="2" & pr$x==2], 0.1125475218393)

  shape <- fs$res$est[fs$res$component==1 & fs$res$terms=="shape"]
  rate <- fs$res$est[fs$res$component==1 & fs$res$terms=="rate"]
  beta <- fs$res$est[fs$res$component==1 & fs$res$terms=="x"]
  rate1 <- exp(log(rate) + beta*1)
  rate2 <- exp(log(rate) + beta*2)
  dens <- pdf_flexsurvmix(fs, newdata=nd, t=1:4)

  expect_equal(dens$dens[dens$event==1 & dens$x==1 & dens$t==1],
               dgamma(1, shape, rate1))
  expect_equal(dens$dens[dens$event==1 & dens$x==2 & dens$t==1],
               dgamma(1, shape, rate2))
  expect_equal(dens$dens[dens$event==1 & dens$x==2 & dens$t==4],
               dgamma(4, shape, rate2))
})

test_that("Aalen Johansen estimates",{
  expect_error(ajfit_flexsurvmix(fs), "Nonparametric estimation not supported with non-factor")
})


