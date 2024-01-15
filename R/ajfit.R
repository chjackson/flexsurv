##' Check the fit of Markov flexible parametric multi-state models against
##' nonparametric estimates.
##'
##' Computes both parametric and comparable Aalen-Johansen nonparametric
##' estimates from a flexible parametric multi-state model, and returns them
##' together in a tidy data frame.  Only models with no covariates, or only
##' factor covariates, are supported.  If there are factor covariates, then the
##' nonparametric estimates are computed for subgroups defined by combinations
##' of the covariates.  Another restriction of this function is that all
##' transitions must have the same covariates on them.
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
##'
##'
##' @return Tidy data frame containing both Aalen-Johansen and parametric
##'   estimates of state occupancy over time, and by subgroup if subgroups are
##'   included.
##'
##' @export
ajfit_fmsm <- function(x, maxt=NULL, newdata=NULL){
  dat <- x[[1]]$data$m
  covnames <- attr(dat,"covnames")
  faccovs <- sapply(dat[,covnames], is.factor)
  if (!all(faccovs)) 
    stop("Nonparametric estimation not supported with non-factor covariates")
  covs <- lapply(x, function(x)attr(x$data$m, "covnames"))
  if (length(covs)>1){
    for (i in 2:length(covs)){
      if(!identical(covs[[i]], covs[[1]]))
        stop("Not currently supported with different covariates on different transitions")
    }
  }
  if (is.null(newdata)) 
    newdata <- do.call(expand.grid, lapply(dat[,covnames,drop=FALSE],  levels))
  else if (is.list(newdata)) newdata <- as.data.frame(newdata)
  else stop("`newdata` should be a data frame")
  
  nmods <- length(x)
  datlist <- vector(nmods, mode="list")
  for (j in 1:nmods){
    datlist[[j]] <-  x[[j]]$data$m
    names(datlist[[j]])[1] <- "(response)"
    datlist[[j]]$trans <- j
  }
  dat <- do.call(dplyr::bind_rows, datlist)
  dat$trans <- factor(dat$trans, labels=attr(x,"names"))
  ## remove interval censored
  dat <- dat[dat$`(response)`[,"status"] != 3,]
  dat$time <- dat[,1][,1]
  dat$status <- dat[,1][,"status"]
  
  ncovvals <- nrow(newdata)
  if (ncovvals==0) ncovvals <- 1
  pt <- vector(ncovvals, mode="list")
  for (i in 1:ncovvals) { 
    datsub <- dat
    for (j in seq_along(covnames)) 
      datsub <- datsub[datsub[,covnames[j]] == newdata[i,covnames[j]],]
    cf <- coxph(Surv(time, status) ~ strata(trans), data=datsub)
    ms <-  msfit(cf, trans=attr(x,"trans"))
    pt[[i]] <- probtrans(ms, predt=0, variance=FALSE)[[1]]
    covvals <- newdata[i,,drop=FALSE][rep(1,nrow(pt[[i]])),,drop=FALSE]
    pt[[i]] <- cbind(pt[[i]], covvals)
  } 
  pt <- do.call(dplyr::bind_rows, pt)
  nstates <- length(attr(x,"statenames"))
  names(pt)[names(pt) %in% paste0("pstate",1:nstates)] <- attr(x,"statenames")
  pt$model <- "Aalen-Johansen"
  
  ## Estimates from parametric competing risks model
  pmat <- vector(ncovvals, mode="list")
  if (is.null(maxt)) 
    maxt <- max(pt$time)
  times <- seq(0, maxt, length.out=100)
  for (i in 1:ncovvals){  # note newdata doesn't support multiple covs 
    pmat[[i]] <- pmatrix.fs(x, newdata=newdata[i,,drop=FALSE], start=1, t=times, tidy=TRUE, ci=FALSE)
    covvals <- newdata[i,,drop=FALSE][rep(1,nrow(pmat[[i]])),,drop=FALSE]
    pmat[[i]] <- cbind(pmat[[i]], covvals)
  } 
  pmat <- do.call(dplyr::bind_rows, pmat)
  pmat <- pmat[pmat$start==1,,drop=FALSE]
  pmat$start <- NULL
  pmat$model <- "Parametric"
  cols <- c("time",covnames,"model",attr(x,"statenames"))
  res <- rbind(pt[,cols], pmat[,cols])
  res <- tidyr::pivot_longer(res, cols=tidyselect::all_of(attr(x,"statenames")), 
                             names_to="state", values_to="val")
  res
}

