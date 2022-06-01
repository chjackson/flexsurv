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


