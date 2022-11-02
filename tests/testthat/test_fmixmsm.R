test_that("fmixmsm state pathways: basic competing risks",{
  tbl_data <- data.frame(
    to = factor(1:3, levels = 1:3, labels = c("response", "progression", "death")),
    dt = rep(1, 3), status = rep(1, 3))
  fit1 <- flexsurvmix(Surv(dt, status) ~ 1, event = tbl_data$to, data = tbl_data,
                      dists = c("response" = "exponential", "progression" = "exponential", "death" = "exponential"))
  
  fit <- fmixmsm("stable" = fit1)
  expect_equivalent(attr(fit, "pathways"),
                    list(c("stable", "response"),
                         c("stable", "progression"),
                         c("stable", "death")))
})

test_that("fmixmsm state pathways: multiple routes to same absorbing state",{
  
  tbl_data <- data.frame(
    to = factor(1:3, levels = 1:3, labels = c("response", "progression", "death")),
    dt = rep(1, 3), status = rep(1, 3))
  fit1 <- flexsurvmix(Surv(dt, status) ~ 1, event = tbl_data$to, data = tbl_data,
                      dists = c("response" = "exponential", "progression" = "exponential", "death" = "exponential"))
  
  tbl_data <- data.frame(
    to = factor(1:2, levels = 1:2, labels = c("progression", "death")),
    dt = rep(1, 2), status = rep(1, 2))
  fit2 <- flexsurvmix(Surv(dt, status) ~ 1, event = tbl_data$to, data = tbl_data,
                      dists = c("progression" = "exponential", "death" = "exponential"))
  
  tbl_data <- data.frame(
    to = factor(1:1, levels = 1:1, labels = c("death")),
    dt = rep(1, 1), status = rep(1, 1))
  fit3 <- flexsurvmix( Surv(dt, status) ~ 1, event = tbl_data$to,  data = tbl_data,
                       dists = c("death" = "exponential"))
  
  fit <- fmixmsm("stable" = fit1,  "response" = fit2,  "progression" = fit3)
  
  expect_equivalent(attr(fit, "pathways"),
                    list(c("stable", "response", "progression", "death"),
                         c("stable", "response", "death"),
                         c("stable", "progression", "death"),
                         c("stable", "death")))
  expect_false(attr(fit, "cycle"))
})

test_that("fmixmsm state pathways: cycles",{
  tbl_data <- data.frame(
    to = factor(1:3, levels = 1:3, labels = c("response", "progression", "death")),
    dt = rep(1, 3), status = rep(1, 3))
  fit1 <- flexsurvmix(Surv(dt, status) ~ 1, event = tbl_data$to, data = tbl_data,
                      dists = c("response" = "exponential", "progression" = "exponential", "death" = "exponential"))
  
  tbl_data <- data.frame(
    to = factor(1:2, levels = 1:2, labels = c("progression", "stable")),
    dt = rep(1, 2), status = rep(1, 2))
  fit2 <- flexsurvmix(Surv(dt, status) ~ 1, event = tbl_data$to, data = tbl_data,
                      dists = c("progression" = "exponential", "stable" = "exponential"))

  fit <- fmixmsm("stable" = fit1,  "response" = fit2)
  expect_true(attr(fit, "cycle"))
})
