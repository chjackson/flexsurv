cuts <- c(1, 3)
haz <- c(0.2, 0.5, 1)        # true hazards in [0,1], (1,3], (3,Inf)
x <- c(0, 0.5, 0.999, 1.001, 2, 3.5, 5, 20)

## Parameter list (rate, hazard ratios) giving hazards `haz` for a baseline
pexp_args <- function(haz, baseline)
  c(list(haz[baseline]), as.list(haz[-baseline] / haz[baseline]))

## Call one of the dfns with the parameters giving hazards `haz`
call_fn <- function(pe, fn, x, haz, baseline = 1, ...)
  do.call(pe$dfns[[fn]], c(list(x), pexp_args(haz, baseline), list(...)))

## Reference values computed directly
H_true <- function(x, haz = c(0.2, 0.5, 1), cuts = c(1, 3)){
  lower <- c(0, cuts); upper <- c(cuts, Inf)
  sapply(x, function(xi) sum(haz * pmax(0, pmin(xi, upper) - lower)))
}

test_that("parameter names depend on the baseline period", {
  expect_equal(custom_pexp(cuts)$dist$pars, c("rate", "hr2", "hr3"))
  expect_equal(custom_pexp(cuts, baseline = 2)$dist$pars, c("rate", "hr1", "hr3"))
  expect_equal(custom_pexp(cuts, baseline = 3)$dist$pars, c("rate", "hr1", "hr2"))
  expect_equal(custom_pexp(numeric(0))$dist$pars, "rate")
})

test_that("dist list has the structure flexsurvreg expects", {
  dist <- custom_pexp(cuts)$dist
  expect_named(dist, c("name", "pars", "location", "transforms",
                       "inv.transforms", "inits"))
  expect_equal(dist$location, "rate")
  expect_length(dist$transforms, 3)
  expect_length(dist$inv.transforms, 3)
  expect_equal(dist$inits(c(1, 2, 3)), c(0.5, 1, 1))
})

test_that("all functions have the parameters as named arguments", {
  pe <- custom_pexp(cuts, baseline = 2)
  for (fn in names(pe$dfns))
    expect_true(all(pe$dist$pars %in% names(formals(pe$dfns[[fn]]))), info = fn)
})

test_that("invalid arguments give errors", {
  expect_error(custom_pexp(c(-1, 2)), "positive")
  expect_error(custom_pexp(c(0, 2)), "positive")
  expect_error(custom_pexp(c(1, Inf)), "positive")
  expect_error(custom_pexp(cuts, baseline = 0), "baseline")
  expect_error(custom_pexp(cuts, baseline = 4), "baseline")
  expect_error(custom_pexp(cuts, baseline = 1:2), "baseline")
})

test_that("change points are sorted and deduplicated", {
  pe1 <- custom_pexp(c(1, 3))
  pe2 <- custom_pexp(c(3, 1, 3))
  expect_equal(call_fn(pe1, "p", x, haz), call_fn(pe2, "p", x, haz))
})

test_that("cumulative hazard, survival and density are correct, for every baseline", {
  for (b in 1:3) {
    pe <- custom_pexp(cuts, baseline = b)
    expect_equal(call_fn(pe, "H", x, haz, b), H_true(x), info = b)
    expect_equal(call_fn(pe, "p", x, haz, b), 1 - exp(-H_true(x)), info = b)
    expect_equal(call_fn(pe, "p", x, haz, b, lower.tail = FALSE),
                 exp(-H_true(x)), info = b)
    expect_equal(call_fn(pe, "h", x, haz, b), haz[c(1, 1, 1, 2, 2, 3, 3, 3)],
                 info = b)
    expect_equal(call_fn(pe, "d", x, haz, b),
                 haz[c(1, 1, 1, 2, 2, 3, 3, 3)] * exp(-H_true(x)), info = b)
  }
})

test_that("results agree with msm::dpexp and ppexp away from change points", {
  skip_if_not_installed("msm")
  pe <- custom_pexp(cuts, baseline = 2)
  expect_equal(call_fn(pe, "d", x, haz, 2), msm::dpexp(x, haz, c(0, cuts)))
  expect_equal(call_fn(pe, "p", x, haz, 2), msm::ppexp(x, haz, c(0, cuts)))
})

test_that("intervals are right-closed at the change points", {
  pe <- custom_pexp(cuts)
  expect_equal(call_fn(pe, "h", c(0, 1, 3), haz), c(0.2, 0.2, 0.5))
})

test_that("density is the derivative of the CDF and integrates to 1", {
  pe <- custom_pexp(cuts)
  eps <- 1e-6
  xx <- c(0.5, 2, 4)
  num_d <- (call_fn(pe, "p", xx + eps, haz) - call_fn(pe, "p", xx - eps, haz)) /
    (2 * eps)
  expect_equal(call_fn(pe, "d", xx, haz), num_d, tolerance = 1e-6)
  total <- integrate(function(t) call_fn(pe, "d", t, haz), 0, Inf)$value
  expect_equal(total, 1, tolerance = 1e-6)
})

test_that("log and log.p options are consistent", {
  pe <- custom_pexp(cuts)
  xx <- x[x > 0]
  expect_equal(call_fn(pe, "d", xx, haz, log = TRUE), log(call_fn(pe, "d", xx, haz)))
  expect_equal(call_fn(pe, "h", xx, haz, log = TRUE), log(call_fn(pe, "h", xx, haz)))
  expect_equal(call_fn(pe, "H", xx, haz, log = TRUE), log(call_fn(pe, "H", xx, haz)))
  expect_equal(call_fn(pe, "p", xx, haz, log.p = TRUE),
               log(call_fn(pe, "p", xx, haz)))
  expect_equal(call_fn(pe, "p", xx, haz, lower.tail = FALSE, log.p = TRUE),
               -H_true(xx))
})

test_that("negative and infinite times are handled", {
  pe <- custom_pexp(cuts)
  expect_equal(call_fn(pe, "d", c(-1, Inf), haz), c(0, 0))
  expect_equal(call_fn(pe, "p", c(-1, Inf), haz), c(0, 1))
  expect_equal(call_fn(pe, "H", c(-1, Inf), haz), c(0, Inf))
})

test_that("quantile function inverts the CDF, for every baseline", {
  pp <- c(0, 0.01, 0.18, 0.3, 0.6, 0.9, 0.99)
  for (b in 1:3) {
    pe <- custom_pexp(cuts, baseline = b)
    qq <- call_fn(pe, "q", pp, haz, b)
    expect_equal(call_fn(pe, "p", qq, haz, b), pp, info = b)
    expect_equal(call_fn(pe, "q", 1 - pp, haz, b, lower.tail = FALSE), qq,
                 info = b)
    expect_equal(call_fn(pe, "q", log(pp[-1]), haz, b, log.p = TRUE), qq[-1],
                 info = b)
  }
  pe <- custom_pexp(cuts)
  expect_equal(call_fn(pe, "q", 1, haz), Inf)
  expect_true(all(is.nan(suppressWarnings(call_fn(pe, "q", c(-0.1, 1.1), haz)))))
})

test_that("parameters are vectorised and recycled against x", {
  pe <- custom_pexp(cuts)
  rates <- c(0.1, 0.2, 0.4)
  vec <- pe$dfns$p(2, rate = rates, hr2 = 2, hr3 = 3)
  one <- sapply(rates, function(r) pe$dfns$p(2, rate = r, hr2 = 2, hr3 = 3))
  expect_equal(vec, one)
  expect_length(pe$dfns$d(c(1, 2, 3, 4), rate = c(0.1, 0.2), hr2 = 2, hr3 = 3), 4)
})

test_that("random generation follows the distribution", {
  pe <- custom_pexp(cuts, baseline = 2)
  set.seed(1)
  sim <- call_fn(pe, "r", 10000, haz, 2)
  expect_length(sim, 10000)
  expect_true(all(sim >= 0))
  ks <- ks.test(sim, function(q) call_fn(pe, "p", q, haz, 2))
  expect_gt(ks$p.value, 0.01)
})

test_that("with no change points the distribution is exponential", {
  pe <- custom_pexp(numeric(0))
  expect_equal(pe$dfns$d(x, rate = 0.3), dexp(x, 0.3))
  expect_equal(pe$dfns$p(x, rate = 0.3), pexp(x, 0.3))
  expect_equal(pe$dfns$q(c(0.1, 0.5, 0.9), rate = 0.3),
               qexp(c(0.1, 0.5, 0.9), 0.3))
})

## Model fitting --------------------------------------------------------------

#skip_if_not_installed("flexsurv")
#library(flexsurv)

bc_split <- survival::survSplit(Surv(recyrs, censrec) ~ ., data = bc,
                                cut = c(2, 4), episode = "int")

fit_pexp <- function(formula, baseline = 1){
  pe <- custom_pexp(c(2, 4), baseline = baseline)
  suppressMessages(flexsurvreg(formula, data = bc, dist = pe$dist,
                               dfns = pe$dfns))
}

test_that("without covariates, fitted hazards equal events / exposure", {
  mle <- with(bc_split, tapply(censrec, int, sum) /
                        tapply(recyrs - tstart, int, sum))
  est <- fit_pexp(Surv(recyrs, censrec) ~ 1)$res[, "est"]
  expect_equal(unname(est["rate"] * c(1, est["hr2"], est["hr3"])),
               unname(as.numeric(mle)), tolerance = 1e-4)
})

test_that("proportional hazards model agrees with Poisson GLM, for every baseline", {
  for (b in 1:3) {
    fit <- fit_pexp(Surv(recyrs, censrec) ~ group, baseline = b)
    bc_split$intf <- relevel(factor(bc_split$int), ref = as.character(b))
    glmfit <- glm(censrec ~ intf + group + offset(log(recyrs - tstart)),
                  family = poisson, data = bc_split)
    expect_equal(unname(fit$res.t[, "est"]), unname(coef(glmfit)),
                 tolerance = 1e-3, info = b)
    expect_equal(fit$loglik, -829.5992, tolerance = 1e-6, info = b)
  }
})

test_that("summary() quantiles match the quantile function", {
  fit <- fit_pexp(Surv(recyrs, censrec) ~ 1)
  est <- as.list(fit$res[, "est"])
  med <- summary(fit, type = "median", ci = FALSE)[[1]]$est
  expect_equal(med, do.call(fit$dfns$q, c(list(0.5), est)))
})
