##' Simulate censored time-to-event data from a fitted flexsurvreg model
##'
##' @param object Object returned by \code{\link{flexsurvreg}}.
##' 
##' @param nsim Number of simulations per row in \code{newdata}.
##' 
##' @param seed Random number seed. This is returned with the result of this
##'   function, as described in \code{\link{simulate}} for the \code{lm} method.
##' 
##' @param newdata Data frame defining alternative sets of covariate values to simulate with.
##' If omitted, this defaults to the data originally used to fit the model. 
##'
##' @param start Delayed entry (left-truncation) time.  If omitted, no delayed entry is assumed.
##'
##' @param censtime A right-censoring time, or vector of times matching the rows
##'   of \code{newdata}.  If \code{NULL} (the default) then uncensored times to events
##' are simulated. 
##'
##' @param tidy If \code{TRUE} then a "tidy" or "long"-format data frame is
##'   returned, with rows defined by combinations of covariates and simulation
##'   replicates.  The simulation replicate is indicated in the column named \code{i}. 
##'
##'   If \code{FALSE}, then a data frame is returned with one row per set of
##'   covariate values, and different columns for different simulation
##'   replicates.  This is the traditional format for `simulate` methods in base
##'   R.
##'
##'   In either case, the simulated time and indicator for whether the time is
##'   an event time (rather than a time of right-censoring) are returned in
##'   different columns.
##'   
##' @param ... Other arguments (not currently used).
##' 
##' @return A data frame, with format determined by whether \code{tidy} was specified.
##' 
##' @examples
##' fit <- flexsurvreg(formula = Surv(futime, fustat) ~ rx, data = ovarian, dist="weibull")
##' fit2 <- flexsurvspline(formula = Surv(futime, fustat) ~ rx, data = ovarian, k=3)
##' nd = data.frame(rx=1:2)
##' simulate(fit, seed=1002, newdata=nd)
##' simulate(fit, seed=1002, newdata=nd, start=500)
##' simulate(fit2, nsim=3, seed=1002, newdata=nd)
##' simulate(fit2, nsim=3, seed=1002, newdata=nd, start=c(500,1000))
##'
##' @export
simulate.flexsurvreg <- function(object, nsim=1, seed=NULL,
                                 newdata=NULL, start=NULL, censtime=NULL, tidy=FALSE,...) {
  if (is.null(newdata))
    newdata <- model.frame(object)
  else if (!is.data.frame(newdata))
    stop("`newdata` should be a data frame")
  if (!is.null(seed)) {
    set.seed(seed)
    attr(seed, "kind") <- as.list(RNGkind())
    seed_attr <- seed 
  } else seed_attr <- .Random.seed
  nd <- nrow(newdata)
  if (!is.numeric(nsim) || (nsim < 1)) stop("`nsim` should be a number >= 1")
  n <- nd*nsim
  if (!is.null(start)) {
    if(!(length(start) %in% c(1, nd)))
      stop(sprintf("`start` of length %s, should be of length 1 or %s = nrow(newdata)"),
           length(start), nd)
    start <- rep(start, length.out=n)
  } else start <- 0
  U <- runif(n, 0, 1)
  newdata <- newdata[rep(1:nd, each=nsim), , drop=FALSE]
  time <- summary(object, newdata=newdata, start=start, 
                 type="quantile", quantiles=U,
                 cross=FALSE, ci=FALSE, se=FALSE, tidy=TRUE)$est
  if (is.null(censtime)) 
    censtime <- rep(Inf, nd)
  if (length(censtime) ==1) 
    censtime <- rep(censtime, nd)
  if (length(censtime)!=nd)
    stop(sprintf("censtime of length %s, should be of length %s = nrow(newdata)"),
         length(censtime), nd)
  censtime <- rep(censtime, each = nsim)
  event <- as.numeric(time <= censtime)
  time <- ifelse(event, time, censtime)
  if (tidy) {
    res <- cbind(newdata, i=rep(1:nsim, nd), time=time, event=event)
    rownames(res) <- NULL
  } else {
    time <- matrix(time, nrow=nd, ncol=nsim, byrow=TRUE)
    event <- matrix(event, nrow=nd, ncol=nsim, byrow=TRUE)
    colnames(time) <- paste0("time_",1:nsim)
    colnames(event) <- paste0("event_",1:nsim)
    res <- as.data.frame(cbind(time, event))
  }
  attr(res, "seed") <- seed_attr
  res
}
