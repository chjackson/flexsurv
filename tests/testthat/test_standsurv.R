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
})

