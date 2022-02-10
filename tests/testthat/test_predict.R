test_that('survival predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")

  # Single time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'survival', t = 500)
  p <- predict(fitw, ovarian, type = 'survival', times = 500)
  expect_equal(s$est, p$.pred)

  # Multiple time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'survival', t = c(500, 1000))
  p <- predict(fitw, ovarian, type = 'survival', times = c(500, 1000))
  expect_equal(s$est, tidyr::unnest(p, .pred)$.pred)
})

test_that('hazard predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")

  # Single time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'hazard', t = 500)
  p <- predict(fitw, ovarian, type = 'hazard', times = 500)
  expect_equal(s$est, p$.pred)

  # Multiple time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'hazard', t = c(500, 1000))
  p <- predict(fitw, ovarian, type = 'hazard', times = c(500, 1000))
  expect_equal(s$est, tidyr::unnest(p, .pred)$.pred)
})

test_that('cumhaz predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")

  # Single time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'cumhaz', t = 500)
  p <- predict(fitw, ovarian, type = 'cumhaz', times = 500)
  expect_equal(s$est, p$.pred)

  # Multiple time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'cumhaz', t = c(500, 1000))
  p <- predict(fitw, ovarian, type = 'cumhaz', times = c(500, 1000))
  expect_equal(s$est, tidyr::unnest(p, .pred)$.pred)
})

test_that('rmst predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")

  # Single time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'rmst', t = 500)
  p <- predict(fitw, ovarian, type = 'rmst', times = 500)
  expect_equal(s$est, p$.pred)

  # Multiple time predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'rmst',
               t = c(500, 1000))
  p <- predict(fitw, ovarian, type = 'rmst', times = c(500, 1000))
  expect_equal(s$est, tidyr::unnest(p, .pred)$.pred)
})

test_that('response/mean predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'mean')
  p <- predict(fitw, ovarian, type = 'response')
  expect_equal(s$est, p$.pred)
})

test_that('quantile predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")

  # Single quantile predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'quantile', t = 500,
               quantiles = c(0.5))
  p <- predict(fitw, ovarian, type = 'quantile', times = 500,
               p = c(0.5))
  expect_equal(s$est, p$.pred)

  # Multiple quantiles predictions
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'quantile',
               quantile = c(0.1, 0.9))
  p <- predict(fitw, ovarian, type = 'quantile', p = c(0.1, 0.9))
  expect_equal(s$est, tidyr::unnest(p, .pred)$.pred)
})

test_that('link predictions', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")
  s <- summary(fitw, ovarian, tidy = TRUE, type = 'link')
  p <- predict(fitw, ovarian, type = 'link')
  expect_equal(s$est, p$.pred)

  p <- predict(fitw, ovarian, type = 'linear')
  expect_equal(s$est, p$.pred)

  p <- predict(fitw, ovarian, type = 'lp')
  expect_equal(s$est, p$.pred)
})

test_that('predictions with missing data', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")
  ovarian_miss <- ovarian
  ovarian_miss$age[[5]] <- NA

  # Single predictions
  p <- predict(fitw, newdata = ovarian_miss)
  expect_equal(nrow(p), nrow(ovarian_miss))

  # Multiple predictions
  p <- predict(fitw, newdata = ovarian_miss,
               type = 'survival', times = c(500, 1000))
  expect_equal(nrow(p), nrow(ovarian_miss))
})

test_that('test order (of age) stays the same', {
  fitw <- flexsurvreg(Surv(futime, fustat) ~ age,
                      data = ovarian, dist = "weibull")
  s <- summary(fitw, newdata = ovarian, type = 'survival', t = 500, tidy = TRUE)
  expect_equal(ovarian$age, s$age)
})
