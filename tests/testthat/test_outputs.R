test_that("normboot.flexsurvreg",{
    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=0))
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
    set.seed(1); b2 <- normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
    expect_equivalent(b1, b2)
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=c(0,50)))
    set.seed(1); b2 <- normboot.flexsurvreg(fite, B=10, X=matrix(c(0,50),nrow=2))
    expect_equivalent(b1, b2)   

    ## return cov effs, not adjusted
    set.seed(1)
    normboot.flexsurvreg(fite, B=5, raw=TRUE)
    set.seed(1)
    normboot.flexsurvreg(fite, B=5, raw=TRUE, transform=TRUE)

    ## no covs
    fite <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
    normboot.flexsurvreg(fite, B=5)
    normboot.flexsurvreg(fite, B=5, transform=TRUE)   
})


test_that("custom function in summary.flexsurvreg",{
    fitw <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="weibull")     

    median.weibull <- function(t, start, shape, scale) { qweibull(0.5, shape=shape, scale=scale) }
    summ <- summary(fitw, newdata=list(age=50), fn=median.weibull, t=1, B=10)
    expect_equal(summ[[1]][1,"est"], 1575.803185910278, tol=1e-04)

    summ <- summary(fitw, newdata=list(age=50), fn=median.weibull, t=c(1,2,3), B=10)
    expect_equal(summ[[1]][1,"est"], summ[[1]][2,"est"])

    mean.weibull <- function(shape, scale=1) { scale * gamma(1 + 1/shape) }   
    median.weibull <- function(t, start, shape, scale) { scale * log(2)^(1/shape) }

})

data(lung)

test_that("newdata in summary.flexsurvreg: dynamic cut, unknown factor level",{
    fl2a <- flexsurvspline(Surv(time, event = status) ~ factor(sex) + cut(age,c(0,56,69,100)), data = lung, k = 2)
    su <- summary(fl2a, newdata = lung, B = 0)
    su1 <- su[[2]][1:5,]
    su <- summary(fl2a, newdata = lung[2,], B = 0)
    su2 <- su[[1]][1:5,]
    expect_equal(su1, su2)
    fl2b <- flexsurvspline(Surv(time, event = status) ~ factor(sex) + cut(age,4), data = lung, k = 2) # should break second summary 
    expect_error(summary(fl2b, newdata = lung[2,], B = 0), "factor .+ has new level")
})

lung$sex <- factor(lung$sex)
fl3 <- flexsurvspline(Surv(time, event = status) ~ sex + age, data = lung, k = 2)

test_that("newdata in summary.flexsurvreg: extra covariates in the list",{
    su1 <- summary(fl3, newdata = lung, B = 0)[[1]][1:5,]
    su2 <- summary(fl3, newdata = lung[1,], B = 0)[[1]][1:5,]
    expect_equal(su1, su2)
})

test_that("newdata in summary.flexsurvreg: missing covariates, factor not supplied as factor",{
    expect_error(summary(fl3, newdata = list(age=10), B = 0), "Value of covariate .+ not supplied")
})

test_that("newdata in summary.flexsurvreg: factor not supplied as factor",{
    lung2 <- lung[1,]; lung2$sex <- as.numeric(1)
    expect_warning(expect_error(summary(fl3, newdata = lung2, B=0), "variable .+ fitted with type"), "not a factor")
    ## numeric doesn't work 
    expect_warning(expect_error(summary(fl3, newdata = list(age=60, sex=1), B=0), "variable .+ fitted with type"), "not a factor")
    ## character works if matches one of the factor levels
    su <- summary(fl3, newdata = list(age=60, sex="1"), B=0) 
    expect_error(summary(fl3, newdata = list(age=60, sex="foo"), B=0), "factor .+ has new level")
})
