##' Check the fit of Markov flexible parametric multi-state models against
##' nonparametric estimates.
##'
##' Computes both parametric and comparable Aalen-Johansen nonparametric
##' estimates from a flexible paramrtric multi-state model, and returns them
##' together in a tidy data frame.  Only models w ith no covariates, or only
##' factor covariates, are supported.  If there are factor covariates, then the
##' nonparametric estimates are computed for subgroups defined by combinations
##' of the covariates.
##'
##' @param x Object returned by \code{\link{fmsm}} representing a flexible
##'   parametric multi-state model.  This must be Markov, rather than
##'   semi-Markov, and no check is performed for this.   Note that all
##'   "competing risks" style models, with just one source state and multiple
##'   destination states, are Markov, so those are fine here.
##'
##' @param maxt Maximum time to compute parametric estimates to.
##'
##' @param newdata Data frame defining the subgroups to consider.  This must
##'   have a column for each covariate in the model.  If omitted, then all
##'   potential subgroups defined by combinations of factor covariates are
##'   included.  Continuous covariates are not supported.
##'   
##' @param ajvar Return standard errors for the Aalen-Johansen estimates. TODO  not implemented yet
##'   
##' @param B Number of bootstrap replicates to use to calculate standard
##'   errors for the parametric estimates.  TODO not implemented yet
##'
##' @return Tidy data frame containing both Aalen-Johansen and parametric
##'   estimates of state occupancy over time, and by subgroup if subgroups are
##'   included.
##'   
ajfit_fmsm <- function(x, maxt=NULL, newdata=NULL, ajvar=FALSE, B=0){
  ## TODO what if different covariates on different transitions? 
  dat <- x[[1]]$data$m
  covnames <- attr(dat,"covnames")
  faccovs <- sapply(dat[,covnames], is.factor)
  if (!all(faccovs)) 
    stop("Nonparametric estimation not supported with non-factor covariates")
  if (is.null(newdata)) 
    newdata <- do.call(expand.grid, lapply(dat[,covnames],  levels))
  
  nmods <- length(x)
  datlist <- vector(nmods, mode="list")
  for (j in 1:nmods){
    datlist[[j]] <-  x[[j]]$data$m
    names(datlist[[j]])[1] <- "(response)"
    datlist[[j]]$trans <- j
  }
  dat <- do.call("rbind", datlist)
  dat$trans <- factor(dat$trans, labels=attr(x,"names"))
  ## remove interval censored
  dat <- dat[dat$`(response)`[,"status"] != 3,]
  dat$time <- dat[,1][,1]
  dat$status <- dat[,1][,"status"]
  
  ncovvals <- nrow(newdata)
  pt <- vector(ncovvals, mode="list")
  for (i in 1:ncovvals) { 
    datsub <- dat
    for (j in seq_along(covnames)) 
      datsub <- datsub[datsub[,covnames[j]] == newdata[i,covnames[j]],]
    cf <- coxph(Surv(time, status) ~ strata(trans), data=datsub)
    ms <-  msfit(cf, trans=attr(x,"trans"))
    pt[[i]] <- probtrans(ms, predt=0, variance=ajvar)[[1]]
    covvals <- newdata[i,][rep(1,nrow(pt[[i]])),]
    pt[[i]] <- cbind(pt[[i]], covvals)
  } 
  pt <- do.call("rbind", pt)
  nstates <- length(attr(x,"statenames"))
  names(pt)[names(pt) %in% paste0("pstate",1:nstates)] <- attr(x,"statenames")
  pt$model <- "Aalen-Johansen"
  
  ## Estimates from parametric competing risks model
  nvals <- nrow(newdata)
  pmat <- vector(nvals, mode="list")
  if (is.null(maxt)) 
    maxt <- max(pt$time)
  times <- seq(0, maxt, length=100)
  for (i in 1:nvals){  # note newdata doesn't support multiple covs 
    if (B>0) ci <- TRUE else ci <- FALSE
    pmat[[i]] <- pmatrix.fs(x, newdata=newdata[i,], start=1, t=times, tidy=TRUE, ci=ci, B=B)
    covvals <- newdata[i,][rep(1,nrow(pmat[[i]])),]
    browser()
    pmat[[i]] <- cbind(pmat[[i]], covvals)
  } 
  pmat <- do.call("rbind", pmat)
  pmat <- pmat[pmat$start==1,,drop=FALSE]
  pmat$start <- NULL
  pmat$model <- "Parametric"
  cols <- c("time",covnames,"model",attr(x,"statenames"))
  res <- rbind(pt[,cols], pmat[,cols])
  res <- tidyr::pivot_longer(res, cols=tidyselect::all_of(attr(x,"statenames")), names_to="state", values_to="prob")
  res
}

