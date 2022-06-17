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


test_that("newdata in summary.flexsurvreg: are missing values passed through or dropped",{
    luna <- lung[1:5,]
    luna$age[1] <- NA
    summ <- summary(fl3, newdata = luna, B = 0, t=100, tidy=TRUE)
    expect_true(is.na(summ$est[1]))
    summ <- summary(fl3, newdata = luna, B = 0, t=100, tidy=TRUE, na.action=na.omit)
    expect_true(!is.na(summ$est[1]))
    
    fitw <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist = "weibull")
    ovarian_miss <- ovarian[1:2,]
    ovarian_miss$age[[1]] <- NA
    summ <- summary(fitw, ovarian_miss, type = 'rmst', t = 500, tidy=TRUE)
    expect_true(is.na(summ$est[1]))
    expect_true(!is.na(summ$est[2]))
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

test_that("summary.flexsurvreg tidy output",{
    lung$sex2 <- as.factor(lung$sex)
    lung$agecat <- ifelse(lung$age<65,"<65",">=65")
    head(lung,20)
    Model.1 <- flexsurvreg(Surv(time, status) ~sex2+agecat ,data=lung, dist="weibull")
    Extrapolation.Data <- model.frame(Model.1)[,c(-1,-dim(model.frame(Model.1))[2])]
    Unique.counts <- as.data.frame(table(Extrapolation.Data))

    ## bug reported by Owain Saunders
    st <- summary.flexsurvreg(Model.1, newdata=Unique.counts[,-3],t= c(4,5),tidy=TRUE, ci=FALSE)
    snt <- summary.flexsurvreg(Model.1, newdata=Unique.counts[,-3],t= c(4,5), ci=FALSE)
    expect_equal(st[st$time==5 & st$sex2==2 & st$agecat=="<65", "est"],
                 snt[["sex2=2,agecat=<65"]][2,2])
    res <- Model.1$res[,"est"]
    res.t <- Model.1$res.t[,"est"]
    expect_equal(pweibull(5, shape=res["shape"], scale=exp(res.t["scale"]+res.t["sex22"]), lower.tail=FALSE),
                 snt[["sex2=2,agecat=<65"]][2,2])

    ## no covariates
    Model.nc <- flexsurvreg(Surv(time, status) ~1 ,data=lung, dist="weibull")
    expect_equivalent(summary.flexsurvreg(Model.nc, t= c(4,5),tidy=TRUE, ci=FALSE),
                      summary.flexsurvreg(Model.nc, t= c(4,5),tidy=FALSE, ci=FALSE)[[1]])

    ## covariates but no newdata - covariate column should be included
    st <- summary.flexsurvreg(Model.1, tidy=TRUE, ci=FALSE)
    expect_equal(st[1,"est"], 0.99604726078272, tol=1e-06)
    expect_equivalent(st[1,"agecat"], ">=65")
})

test_that("summary.flexsurvreg untidy output back compatibility",{
  nd <- lung[c(1,1,2),]
  s1 <- summary(fl3, newdata = nd, B=0, tidy=TRUE, t=c(5,10))
  s2 <- summary(fl3, newdata = nd, B=0, tidy=FALSE, t=c(5,10))
  expect_equal(s1[3,"est"], s2[[2]][1,"est"])  
})

test_that("hazard ratio",{
    t <- c(100,200,300)
    haz2 <- summary(fl3, type="hazard", t=t, newdata=list(age=60, sex="2"), ci=FALSE, tidy=TRUE)
    haz1 <- summary(fl3, type="hazard", t=t, newdata=list(age=60, sex="1"), ci=FALSE, tidy=TRUE)
    hr <- haz2$est / haz1$est   
    
    hr2 <- hr_flexsurvreg(fl3, newdata=as.data.frame(list(age=60, sex=c("1","2"))), t=t, ci=FALSE)
    expect_equal(hr, hr2$est)

    hr2 <- hr_flexsurvreg(fl3, newdata=as.data.frame(list(age=60, sex=c("1","2"))), 
                          t=t, ci=TRUE, B=10)
    expect_is(hr2$lcl, "numeric")
    
    ## a non-proportional hazards model
    fl4 <- flexsurvspline(Surv(time, event = status) ~ sex + age, anc=list(gamma1=~sex), 
                          data = lung, k = 2)
    hr4 <- hr_flexsurvreg(fl4, newdata=as.data.frame(list(age=60, sex=c("1","2"))), 
                          t=t, ci=TRUE, B=10)
    expect_false(hr4$est[1] == hr4$est[2])

    expect_error(hr_flexsurvreg(fl4, t=t), "`newdata` must be specified")

    ## default newdata: factor
    fitw <- flexsurvreg(Surv(futime, fustat) ~ factor(rx), data = ovarian, dist="weibull")
    expect_equal(hr_flexsurvreg(fitw, t=t)$est,
                 hr_flexsurvreg(fitw, t=t, newdata=list(rx=c("1","2")))$est)
    ## default newdata: numeric
    ovarian$rxbin <- ovarian$rx - 1 
    fitw <- flexsurvreg(Surv(futime, fustat) ~ rxbin, data = ovarian, dist="weibull")
    expect_equal(hr_flexsurvreg(fitw, t=t)$est,
                 hr_flexsurvreg(fitw, t=t, newdata=list(rxbin=c(0,1)))$est)
})
