test_that('survival predictions', {
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ group, 
                      data=bc, dist="weibull")
  
  # Single time predictions
  ss <- standsurv(fitw, at=list(list(group="Good"), 
                                            list(group="Medium"), 
                                            list(group="Poor")),
                              t = 4)
  s <- summary(fitw, tidy = TRUE, type = 'survival', t = 4, ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est)

  # Multiple time predictions
  ss <- standsurv(fitw, at=list(list(group="Good"), 
                                            list(group="Medium"), 
                                            list(group="Poor")),
                              t = c(4,5))
  s <- summary(fitw, tidy = TRUE, type = 'survival', t = c(4,5), ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est)
  
})

test_that('hazard predictions', {
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ group, 
                      data=bc, dist="weibull")
  
  # Single time predictions
  ss <- standsurv(fitw, type="hazard",
                              at=list(list(group="Good"), 
                                      list(group="Medium"), 
                                      list(group="Poor")),
                              t = 4)
  s <- summary(fitw, tidy = TRUE, type = 'hazard', t = 4, ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est)
  
  # Multiple time predictions
  ss <- standsurv(fitw, type="hazard",
                              at=list(list(group="Good"), 
                                      list(group="Medium"), 
                                      list(group="Poor")),
                              t = c(4,5))
  s <- summary(fitw, tidy = TRUE, type = 'hazard', t = c(4,5), ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est)
})

test_that('quantile predictions', {
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ group, 
                      data=bc, dist="weibull")
  
  # Multiple quantile predictions (q=0.1 and 0.5)
  ss <- standsurv(fitw, type="quantile",
                              at=list(list(group="Good"), 
                                      list(group="Medium"), 
                                      list(group="Poor")),
                  quantiles = seq(0.1,0.9, by=0.1))
  s <- summary(fitw, tidy = TRUE, type = 'quantile', quantiles = seq(0.1,0.9, by=0.1), 
               ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est, 
               tolerance = .Machine$double.eps^0.25) # use same tolerance as uniroot
  
  # Test marginal predictions of quantiles correspond to inverse marginal survival
  set.seed(136)
  bc$age <- rnorm(dim(bc)[1], mean = 65 - bc$recyrs, sd = 5)
  fitw2 <- flexsurvreg(Surv(recyrs, censrec) ~ group + age, 
                      data=bc, dist="weibull")
  ss <- standsurv(fitw2, type="quantile",
                  at=list(list(group="Good")),
                  quantiles = seq(0.1, 0.9, by=0.1))
  ss
  # Feed back in the quantiles and calculate marginal survival probabilities
  ss2 <- standsurv(fitw2, type="survival",
                   at=list(list(group="Good")),
                   t=ss$at1)
  expect_equal(1-seq(0.1, 0.9, by=0.1), ss2$at1, 
               tolerance = .Machine$double.eps^0.25) # use same tolerance as uniroot
})

test_that('rmst predictions', {
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ group, 
                      data=bc, dist="weibull")
  
  # Single time predictions
  ss <- standsurv(fitw, type="rmst",
                  at=list(list(group="Good"), 
                          list(group="Medium"), 
                          list(group="Poor")),
                  t = 4)
  s <- summary(fitw, tidy = TRUE, type = 'rmst', t = 4, ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est)
  
  # Multiple time predictions
  ss <- standsurv(fitw, type="rmst",
                  at=list(list(group="Good"), 
                          list(group="Medium"), 
                          list(group="Poor")),
                  t = c(4,5))
  s <- summary(fitw, tidy = TRUE, type = 'rmst', t = c(4,5), ci = FALSE)
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s$est)
})


test_that('marginal_predictions', {
  set.seed(136)
  bc$age <- rnorm(dim(bc)[1], mean = 65 - bc$recyrs, sd = 5)
  
  # Single time marginal predictions for survival
  fitw2 <- flexsurvreg(Surv(recyrs, censrec) ~ group + age, 
                       data=bc, dist="weibull")
  ss <- standsurv(fitw2, at=list(list(group="Good"), 
                                             list(group="Medium"), 
                                             list(group="Poor")),
                              t = 4)
  s.good <- summary(fitw2, bc %>% mutate(group="Good"), tidy = TRUE, 
                    type = 'survival', t = 4, ci = FALSE)
  s.medium <- summary(fitw2, bc %>% mutate(group="Medium"), tidy = TRUE, 
                      type = 'survival', t = 4, ci = FALSE)
  s.poor <- summary(fitw2, bc %>% mutate(group="Poor"), tidy = TRUE, 
                    type = 'survival', t = 4, ci = FALSE)
  s.marginal <- c(mean(s.good$est), mean(s.medium$est), mean(s.poor$est))
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s.marginal)
  
  # Single time marginal predictions for hazard
  ss <- standsurv(fitw2, type="hazard",
                              at=list(list(group="Good"), 
                                      list(group="Medium"), 
                                      list(group="Poor")),
                              t = 4)
  h.good <- summary(fitw2, bc %>% mutate(group="Good"), tidy = TRUE, 
                    type = 'hazard', t = 4, ci = FALSE)
  h.medium <- summary(fitw2, bc %>% mutate(group="Medium"), tidy = TRUE, 
                      type = 'hazard', t = 4, ci = FALSE)
  h.poor <- summary(fitw2, bc %>% mutate(group="Poor"), tidy = TRUE, 
                    type = 'hazard', t = 4, ci = FALSE)
  s.marginal <- c(weighted.mean(h.good$est, s.good$est), 
                  weighted.mean(h.medium$est, s.medium$est), 
                  weighted.mean(h.poor$est, s.poor$est))
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s.marginal)
  
  # Single time marginal predictions for rmst
  ss <- standsurv(fitw2, type="rmst",
                              at=list(list(group="Good"), 
                                      list(group="Medium"), 
                                      list(group="Poor")),
                              t = 4)
  rmst.good <- summary(fitw2, bc %>% mutate(group="Good"), tidy = TRUE, 
                       type = 'rmst', t = 4, ci = FALSE)
  rmst.medium <- summary(fitw2, bc %>% mutate(group="Medium"), tidy = TRUE, 
                         type = 'rmst', t = 4, ci = FALSE)
  rmst.poor <- summary(fitw2, bc %>% mutate(group="Poor"), tidy = TRUE, 
                       type = 'rmst', t = 4, ci = FALSE)
  s.marginal <- c(mean(rmst.good$est), 
                  mean(rmst.medium$est), 
                  mean(rmst.poor$est))
  expect_equal(c(ss$at1, ss$at2, ss$at3) ,s.marginal)
})


test_that('predictions with missing data', {
  set.seed(136)
  bc$age <- rnorm(dim(bc)[1], mean = 65 - bc$recyrs, sd = 5)
  bc_miss <- bc
  bc_miss$age[[5]] <- NA
  
  fitw2 <- flexsurvreg(Surv(recyrs, censrec) ~ group + age, 
                       data=bc_miss, dist="weibull")
  # Single time predictions
  expect_warning(standsurv(fitw2, type="rmst",
                              at=list(list(group="Good"), 
                                      list(group="Medium"), 
                                      list(group="Poor")),
                              t = 4,
                              newdata = bc_miss))
})


test_that('model with no covariates', {
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ 1, 
                       data=bc, dist="weibull")
  # Single time prediction - checking survival with predict.flexsurvreg
  expect_equal(as.numeric(standsurv(fitw, t = 4, newdata = bc)[1,2]), 
               as.numeric(predict(fitw, type = "survival", times = 4, 
                                  newdata = bc)[1,2]))
})

test_that('all-cause predictions from RS model', {
  # Single time prediction - checking all-cause survival from an RS model
  # with predict.flexsurvreg, where background rates are zero
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ 1, 
                      data=bc, dist="weibull")
  set.seed(136)
  bc$age <- rnorm(dim(bc)[1], mean = 65 - bc$recyrs, sd = 5)
  bc$agedays <- floor(bc$age * 365.25)
  ## Create some random diagnosis dates centred on 01/01/2010 with SD=1 year
  bc$diag <- as.Date(floor(rnorm(dim(bc)[1], mean = as.Date("01/01/2010", "%d/%m/%Y"), sd=365)), origin="1970-01-01")
  ## Create sex (assume all are female)
  bc$sex <- factor("female")
  ## Assume background hazard of zero
  bc$bhazard <- 0
  ## ratetable (used by standsurv)
  new.ratetable <- survexp.us 
  new.ratetable[!is.na(new.ratetable)] <- 0
  fitw2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, 
                       data=bc, dist="weibull", bhazard=bhazard)
  pred1 <- as.numeric(predict(fitw, type = "survival", times = 4, 
                              newdata = bc)[1,2])
  pred2 <- as.numeric(standsurv(fitw2, t=4, newdata=bc, type="survival",
            rmap=list(sex = sex,
                      year = diag,
                      age = agedays
            ),
            ratetable = new.ratetable,
            scale.ratetable = 365.25)[1,2])
  expect_equal(pred1, pred2)
  
  # all-cause survival quantile
  # just use one row of data
  pred1 <- summary(fitw, type = "quantile", quantiles = seq(0.1, 0.9, 0.1), 
                              newdata = bc[1,], ci=F, tidy=T)
  pred2 <- standsurv(fitw2, quantiles = seq(0.1, 0.9, 0.1), 
                                newdata=bc[1,], type="quantile",
                                rmap=list(sex = sex,
                                          year = diag,
                                          age = agedays
                                ),
                                ratetable = new.ratetable,
                                scale.ratetable = 365.25)
  expect_equal(pred1$est, pred2$at1, tolerance = .Machine$double.eps^0.25) # use same tolerance as uniroot)
  
})

test_that("standsurv deltamethod and bootstrap SEs",{
  fitw <- flexsurvreg(Surv(recyrs, censrec) ~ group, 
                      data=bc, dist="weibull")
  ss <- standsurv(fitw, at=list(list(group="Good"), list(group="Medium"), list(group="Poor")),
                  t = 4, se=TRUE, boot=FALSE)
  expect_true(is.numeric(ss$at1_se))
  ssb <- standsurv(fitw, at=list(list(group="Good"), list(group="Medium"), list(group="Poor")),
                  t = 4, se=TRUE, boot=TRUE, B=5)
  expect_true(ssb$at1_se != ss$at1_se)
})

if (interactive() || covr::in_covr()){
  test_that("standsurv plot",{
    expect_error({
                 
    ## Use bc dataset, with an age variable appended
    ## mean age is higher in those with smaller observed survival times 
    newbc <- bc
    newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE), 
                       sd = 5)

    ## Fit a Weibull flexsurv model with group and age as covariates
    weib_age <- flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc,
                            dist="weibull")
    ## Calculate standardized survival and the difference in standardized survival
    ## for the three levels of group across a grid of survival times
    standsurv_weib_age <- standsurv(weib_age,
                                    at = list(list(group="Good"),
                                              list(group="Medium"),
                                              list(group="Poor")),
                                    t=seq(0,7, length=100),
                                    contrast = "difference", ci=TRUE,
                                    boot = TRUE, B=10, seed=123)
    plot(standsurv_weib_age)
    plot(standsurv_weib_age) + ggplot2::theme_bw() + ggplot2::ylab("Survival") +
      ggplot2::xlab("Time (years)") + 
      ggplot2::guides(color=ggplot2::guide_legend(title="Prognosis"),
                      fill=ggplot2::guide_legend(title="Prognosis"))
    plot(standsurv_weib_age, contrast=TRUE, ci=TRUE) + 
      ggplot2::ylab("Difference in survival") 
    }, NA)
  })
}

test_that("examples from help(standsurv) run",{
  skip_on_cran()
  expect_error({
    ## mean age is higher in those with smaller observed survival times 
    newbc <- bc
    set.seed(1)
    newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE),
                       sd = 5)

    ## Fit a Weibull flexsurv model with group and age as covariates
    weib_age <- flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc, 
                            dist="weibull")

    ## Calculate standardized survival and the difference in standardized survival
    ## for the three levels of group across a grid of survival times                        
    standsurv_weib_age <- standsurv(weib_age, 
                                    at = list(list(group="Good"), 
                                              list(group="Medium"), 
                                              list(group="Poor")), 
                                    t=seq(0,7, length.out=100),
                                    contrast = "difference", ci=FALSE)
    standsurv_weib_age

    ## Calculate hazard of standardized survival and the marginal hazard ratio
    ## for the three levels of group across a grid of survival times
    ## 10 bootstraps for confidence intervals (this should be larger)

    haz_standsurv_weib_age <- standsurv(weib_age, 
                                        at = list(list(group="Good"), 
                                                  list(group="Medium"), 
                                                  list(group="Poor")), 
                                        t=seq(0,7, length.out=100),
                                        type="hazard",
                                        contrast = "ratio", boot = TRUE,
                                        B=10, ci=TRUE)
    haz_standsurv_weib_age

    if (interactive()){
      plot(haz_standsurv_weib_age, ci=TRUE)
      ## Hazard ratio plot shows a decreasing marginal HR 
      ## Whereas the conditional HR is constant (model is a PH model)
      plot(haz_standsurv_weib_age, contrast=TRUE, ci=TRUE)
    }
    ## Calculate standardized survival from a Weibull model together with expected
    ## survival matching to US lifetables

    ## age at diagnosis in days. This is required to match to US ratetable, whose
    ## timescale is measured in days
    newbc$agedays <- floor(newbc$age * 365.25)  
    ## Create some random diagnosis dates centred on 01/01/2010 with SD=1 year
    ## These will be used to match to expected rates in the lifetable
    newbc$diag <- as.Date(floor(rnorm(dim(newbc)[1], 
                                      mean = as.Date("01/01/2010", "%d/%m/%Y"), sd=365)), 
                          origin="1970-01-01")
    ## Create sex (assume all are female)
    newbc$sex <- factor("female")
    standsurv_weib_expected <- standsurv(weib_age, 
                                         at = list(list(group="Good"), 
                                                   list(group="Medium"), 
                                                   list(group="Poor")), 
                                         t=seq(0,7, length.out=100),
                                         rmap=list(sex = sex,
                                                   year = diag,
                                                   age = agedays),
                                         ratetable = survival::survexp.us,
                                         scale.ratetable = 365.25,
                                         newdata = newbc)
    ## Plot marginal survival with expected survival superimposed                                            
    plot(standsurv_weib_expected, expected=TRUE)
  }, NA)
})
