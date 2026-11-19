##' Piecewise exponential distribution for use with flexsurvreg()
##'
##' @param cuts Times of the change points (positive, any order).  These
##'   define \code{K = length(cuts) + 1} periods, \code{[0, cuts[1]], (cuts[1], cuts[2]],
##'   ..., (cuts[K-1], Inf)}.
##'
##' @param baseline Index of the period whose hazard is the baseline
##'   hazard \code{rate}.  The other parameters are hazard ratios
##'   named by period, so that the hazard in period \code{k} is
##'   \code{rate * hrk}.  For example with K = 3 and \code{baseline =
##'   1} (the default) the parameters are \code{rate}, \code{hr2},
##'   \code{hr3}; with \code{baseline = 2} they are \code{rate},
##'   \code{hr1}, \code{hr3}.  \code{rate} is the location parameter,
##'   so covariates on it give a proportional hazards model.
##'
##' @return A list with components
##'   \code{dist}: the custom distribution list, to pass as the \code{dist}
##'   argument of \code{flexsurvreg}, and
##'   \code{dfns}: the list of d, p, h, H, q, r functions, to pass as the
##'   \code{dfns} argument.
##'
##' @examples
##' pe <- custom_pexp(c(365, 730))
##' flexsurvreg(Surv(recyrs*365, censrec) ~ 1, data = bc,
##'             dist = pe$dist, dfns = pe$dfns)
##'
##' @export
custom_pexp <- function(cuts, baseline = 1){
  cuts <- sort(unique(cuts))
  if (any(!is.finite(cuts)) || any(cuts <= 0))
    stop("change points must be positive and finite")
  K <- length(cuts) + 1
  if (length(baseline) != 1 || !(baseline %in% seq_len(K)))
    stop(sprintf("`baseline` must be a single integer between 1 and %d", K))
  pars <- c("rate", if (K > 1) paste0("hr", seq_len(K)[-baseline]))
  lower <- c(0, cuts)
  upper <- c(cuts, Inf)

  ## Recycle x and the parameters (rate, hazard ratios) to a common length n.
  ## Returns x (length n) and R (n x K matrix of hazards in each interval)
  recycle <- function(x, args){
    n <- max(length(x), lengths(args))
    P <- matrix(unlist(lapply(args, rep_len, n)), nrow = n, ncol = K)
    R <- matrix(P[, 1], nrow = n, ncol = K)
    if (K > 1) R[, -baseline] <- P[, 1] * P[, -1]
    list(x = rep_len(x, n), R = R)
  }

  ## Cumulative hazard
  cumhaz <- function(x, R){
    n <- length(x)
    X <- matrix(x, n, K)
    E <- pmax(0, pmin(X, matrix(upper, n, K, byrow = TRUE)) -
                 matrix(lower, n, K, byrow = TRUE))
    E[R == 0 & E == Inf] <- 0
    rowSums(R * E)
  }

  ## Hazard (zero for negative times).  Intervals are right-closed, so an
  ## event exactly at a change point gets the hazard of the earlier interval,
  ## as in survival::survSplit()
  haz <- function(x, R){
    j <- findInterval(x, lower, left.open = TRUE)
    j[x == 0] <- 1
    h <- numeric(length(x))
    pos <- j > 0
    h[pos] <- R[cbind(which(pos), j[pos])]
    h
  }

  ## Each function has formals rate, hr..., as flexsurvreg requires.
  ## mget(pars) collects these from the function's evaluation frame.
  rate_formals <- setNames(rep(list(1), K), pars)
  with_rates <- function(f, before, after){
    formals(f) <- c(before, rate_formals, after)
    f
  }

  d <- with_rates(function(x, log = FALSE){
    a <- recycle(x, mget(pars))
    ld <- log(haz(a$x, a$R)) - cumhaz(a$x, a$R)
    if (log) ld else exp(ld)
  }, alist(x = ), alist(log = FALSE))

  p <- with_rates(function(q, lower.tail = TRUE, log.p = FALSE){
    a <- recycle(q, mget(pars))
    logS <- -cumhaz(a$x, a$R)
    if (lower.tail) {
      if (log.p) log(-expm1(logS)) else -expm1(logS)
    } else {
      if (log.p) logS else exp(logS)
    }
  }, alist(q = ), alist(lower.tail = TRUE, log.p = FALSE))

  h <- with_rates(function(x, log = FALSE){
    a <- recycle(x, mget(pars))
    res <- haz(a$x, a$R)
    if (log) log(res) else res
  }, alist(x = ), alist(log = FALSE))

  H <- with_rates(function(x, log = FALSE){
    a <- recycle(x, mget(pars))
    res <- cumhaz(a$x, a$R)
    if (log) log(res) else res
  }, alist(x = ), alist(log = FALSE))

  q <- with_rates(function(p, lower.tail = TRUE, log.p = FALSE){
    a <- recycle(p, mget(pars))
    p <- a$x
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    Ht <- -log1p(-p)                        # target cumulative hazard
    n <- length(p)
    ## cumulative hazard at the start of each interval (n x K)
    Hlow <- matrix(0, n, K)
    for (k in seq_len(K)[-1])
      Hlow[, k] <- Hlow[, k - 1] + a$R[, k - 1] * (lower[k] - lower[k - 1])
    j <- pmax(1, rowSums(Hlow <= Ht))
    idx <- cbind(seq_len(n), j)
    res <- lower[j] + (Ht - Hlow[idx]) / a$R[idx]
    res[p < 0 | p > 1 | is.na(p)] <- NaN
    res[p == 1] <- Inf
    res
  }, alist(p = ), alist(lower.tail = TRUE, log.p = FALSE))

  r <- with_rates(function(n){
    rates <- mget(pars)
    if (length(n) > 1) n <- length(n)
    do.call(q, c(list(p = runif(n)), rates))
  }, alist(n = ), list())

  dist <- list(
    name = "pexp_custom",
    pars = pars,
    location = "rate",
    transforms = rep(list(log), K),
    inv.transforms = rep(list(exp), K),
    inits = function(t, ...) c(1 / mean(t), rep(1, K - 1))
  )

  list(dist = dist, dfns = list(d = d, p = p, h = h, H = H, q = q, r = r))
}
