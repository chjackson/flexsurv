##' Mean times to events from a flexsurvmix model
##'
##' This returns the mean of each event-specific parametric time-to-event
##' distribution in the mixture model, which is the mean time to event
##' conditionally on that event being the one that happens.
##'
##' @param x Fitted model object returned from \code{\link{flexsurvmix}}.
##'
##' @param newdata Data frame or list of covariate values.   If omitted for a
##'   model with covariates, a default is used, defined by all combinations of
##'   factors if the only covariates in the model are factors, or all covariate
##'   values of zero if there are any non-factor covariates in the model.
##'
##' @param B Number of simulations to use to compute 95\% confidence intervals,
##'   based on the asymptotic multivariate normal distribution of the basic
##'   parameter estimates.  If \code{B=NULL} then intervals are not computed.
##'
##' @return Mean times to next event conditionally on each alternative event,
##'   given the specified covariate values.
##'
##' @export
mean_flexsurvmix <- function(x, newdata=NULL, B=NULL){
  if (is.null(newdata)) newdata <- default_newdata(x)
  if (!is.null(newdata)){  # mean functions aren't vectorised 
    newdata <- as.data.frame(newdata)
    nvals <- nrow(newdata)
    res <- cisumm_flexsurvmix(x, newdata=newdata[1,,drop=FALSE], fnname="mean", fnlist=x$dfns, B=B)
    if (nvals > 1) {
      for (i in 2:nvals){
        res <- rbind(res, cisumm_flexsurvmix(x, newdata=newdata[i,,drop=FALSE], fnname="mean", fnlist=x$dfns, B=B))
      }
    }
  } else {
    res <- cisumm_flexsurvmix(x, newdata=newdata, fnname="mean", fnlist=x$dfns, B=B)
  }
  res
}

##' Restricted mean times to events from a flexsurvmix model
##'
##' This returns the restricted mean of each event-specific parametric time-to-event
##' distribution in the mixture model, which is the mean time to event
##' conditionally on that event being the one that happens, and conditionally
##' on the event time being less than some time horizon \code{tot}.
##'   
##' @param tot Time horizon to compute the restricted mean until.
##'
##' @inheritParams mean_flexsurvmix
##'
##' @return Restricted mean times to next event conditionally on each alternative event,
##'   given the specified covariate values.
##'   
##' @export
rmst_flexsurvmix <- function(x, newdata=NULL, tot=Inf, B=NULL){
  if (is.null(newdata)) newdata <- default_newdata(x)
  if (!is.null(newdata)){  # mean functions aren't vectorised 
    newdata <- as.data.frame(newdata)
    nvals <- nrow(newdata)
    res <- cisumm_flexsurvmix(x, newdata=newdata[1,,drop=FALSE], fnname="rmst", fnarg="t", fnargval=tot, fnlist=x$dfns, B=B)
    if (nvals > 1) {
      for (i in 2:nvals){
        res <- rbind(res, cisumm_flexsurvmix(x, newdata=newdata[i,], fnname="rmst", fnarg="t", fnargval=tot, fnlist=x$dfns, B=B))
      }
    }
  } else {
    newdata <- default_newdata(x, zero=TRUE)
    res <- cisumm_flexsurvmix(x, newdata=newdata, fnname="rmst", fnarg="t", fnargval=tot, fnlist=x$dfns, B=B)
  }
  res
}


##' Quantiles of time-to-event distributions in a flexsurvmix model
##' 
##' This returns the quantiles of each event-specific parametric time-to-event
##' distribution in the mixture model, which describes the time to the event
##' conditionally on that event being the one that happens.
##' 
##' @inheritParams mean_flexsurvmix
##' 
##' @param probs Vector of alternative quantiles, by default \code{c(0.025, 0.95, 0.975)}
##' giving the median and a 95\% interval. 
##' 
##' @export
quantile_flexsurvmix <- function(x, newdata=NULL, B=NULL, probs=c(0.025, 0.5, 0.975)){
  cisumm_flexsurvmix(x, newdata=newdata, fnname="q", fnarg="p", fnargval=probs, fnlist=x$dfns, B=B)
  
}

##' Transition probabilities from a flexsurvmix model
##'
##' These quantities are variously known as transition probabilities, or state
##' occupancy probabilities, or values of the "cumulative incidence" function,
##' or values of the "subdistribution" function. They are the probabilities that
##' an individual has experienced an event of a particular kind by time
##' \code{t}.
##'
##' Note that "cumulative incidence" is a misnomer, as "incidence" typically
##' means a hazard, and the quantities computed here are not cumulative hazards,
##' but probabilities.
##'
##'
##' @inheritParams mean_flexsurvmix
##'
##' @param t Vector of times \code{t} to calculate the probabilities of
##'   transition by.
##'
##' @param startname Name of the state where individuals start. This considers
##'   the model as a multi-state model where people start in this state, and may
##'   transition to one of the competing events.
##'
##' @return A data frame with transition probabilities by time, covariate value
##'   and destination state.
##'
##'
##' @export
p_flexsurvmix <- function(x, newdata=NULL,startname="start", t=1, B=NULL){
  fm_p_fn <- function(base,prob,resdf,x,startname) {
    resdf$val <- base * prob
    resdf$tval <- resdf$pval <- NULL
    resdf  <- resdf %>% 
      tidyr::pivot_wider(names_from="event", values_from="val")
    resdf[[startname]] <- 1 - rowSums(resdf[,x$evnames]) 
    statenames <- c(startname,x$evnames)
    resdf %>% 
      pivot_longer(cols=all_of(statenames), names_to="state", values_to="val")
  }
  res <- cisumm_flexsurvmix(x,   newdata=newdata,  fnname="p", fnarg="q",
                            parclass = c("time", "prob"),
                            fnargval=t,  fnlist=x$dfns, combfn = fm_p_fn, B=B, 
                            startname=startname)
  names(res)[names(res)=="q"] <- "t"
  res
}



##' Probabilities of competing events from a flexsurvmix model
##'
##' @inheritParams mean_flexsurvmix
##'
##' @return A data frame containing the probability that each of the competing
##'   events will occur next, by event and by any covariate values specified in
##'   \code{newdata}.
##'
##' @export
probs_flexsurvmix <- function(x, newdata=NULL, B=NULL){
  if (is.null(newdata)) newdata <- default_newdata(x)
  if (x$K==1) { 
    res <- newdata
    res$event <- x$evnames
    res$val <- 1
    
  } else 
    res <- cisumm_flexsurvmix(x, newdata=newdata, parclass="prob", B=B, fnname=NULL)
  res
}


default_newdata <- function(x, zero=FALSE){
  dat <- x$data$mfcomb
  covnames <- attr(dat, "covnames")
  ncovs <- length(covnames)
  if (ncovs == 0){
    newdata <- NULL
  }
  else { 
    faccovs <- sapply(dat[,covnames], is.factor)
    if (!all(faccovs) & !zero) 
      newdata <- matrix(0, nrow=1, ncol=ncovs, dimnames=list(NULL, covnames))
    else 
      newdata <- do.call(expand.grid, lapply(dat[,covnames],  levels))
    newdata <- as.data.frame(newdata)
  }
  newdata
}

cisumm_flexsurvmix <- function(x, newdata=NULL, 
                               parclass = "time", 
                               fnname, fnarg=NULL, fnargval=NULL, fnlist=NULL, combfn=NULL,
                               B=NULL, ...){
  K <- x$K
  res <- vector(K, mode="list")
  if (is.null(newdata)) newdata <- default_newdata(x)
  if ("time" %in% parclass) { 
    nquants <- max(1, length(fnargval))
    if (!is.null(newdata)) {
      newdata <- as.data.frame(newdata)
      ncovs <- nrow(newdata)
      newdatarep <- newdata[rep(1:ncovs, nquants),,drop=FALSE]
    } else {
      ncovs <- 1 
      newdatarep <- NULL
    }
    fnargvalrep <- rep(fnargval, each=ncovs)
    for (k in 1:x$K){ 
      pars <- get_basepars(x, newdata=newdatarep, event=k)
      parnames <- names(pars)
      if (!is.null(fnarg)) 
        pars[[fnarg]] <- fnargvalrep
      if (!is.null(fnlist))
        fn <- fnlist[[k]][[fnname]]
      res[[k]] <- do.call(fn, pars)
    } 
    tdf <- data.frame(event = rep(x$evnames, each=length(res[[1]])))
    if (!is.null(fnarg)) 
      tdf[[fnarg]] <- rep(fnargvalrep, K)
    if (!is.null(newdata))
      tdf <- cbind(tdf, newdatarep[rep(1:nrow(newdatarep), K),,drop=FALSE])
    tdf$val <- unlist(res)
  } 
  if ("prob" %in% parclass) {   
    if (!is.null(newdata)) newdata <- as.data.frame(newdata)
    pdf <- get_probpars(x, newdata=newdata)
    if (!is.null(newdata)){
      if (length(intersect(names(pdf), names(newdata)) > 0))
        pdf <- dplyr::left_join(pdf, newdata, by=intersect(names(pdf), names(newdata)))
      else pdf <- tidyr::crossing(pdf, newdata)
    } 
  }
  if (setequal(parclass, c("time","prob"))){
    tdf <- tdf %>% dplyr::rename(tval="val") 
    pdf <- pdf %>% dplyr::rename(pval="val") 
    resdf <- tdf %>% dplyr::left_join(pdf, by=c("event", names(newdata)))
    resdf <- combfn(resdf$tval, resdf$pval, resdf, x, ...)
    resdf$tval <- resdf$pval <- NULL
  }  else if (parclass=="time"){ 
    resdf  <- tdf
  }  else if (parclass=="prob") {
    resdf <- pdf
  }
  if (is.numeric(B)) {
    resm <- array(dim=c(B, nrow(resdf)))
    resm[1,] <- as.matrix(resdf$val)
    if (B > 2) {
      for (b in 2:B){
        pars <- list(x = resample_pars(x), newdata=newdata, parclass=parclass, 
                     fnname=fnname, fnarg=fnarg, fnargval=fnargval, fnlist=fnlist, combfn=combfn,
                     B = NULL)
        pars <- c(pars, list(...))
        resm[b,]<- do.call(cisumm_flexsurvmix, pars)$val
      }
    }
    qp <- c(0.025, 0.975)
    ci <- apply(resm, 2, quantile, qp)
    resdf$lower <- ci["2.5%",]
    resdf$upper <- ci["97.5%",]
  }
  resdf
}

resample_pars <- function(x){ 
  cov <- matrix(0, x$npars, x$npars)
  cov[x$optpars, x$optpars] <- x$cov
  newpars <- rmvnorm(1, x$res$est.t[-1], cov)
  x$res$est.t <- c(NA, newpars)
  x$res$est <- inv.transform.res(x, x$dlist) 
  x
}


##' Evaluate baseline time-to-event distribution parameters given covariate values in a flexsurvmix model
##' 
##' @param x Fitted model object
##' 
##' @param newdata Data frame of alternative covariate values 
##' 
##' @param event Event 
##' 
get_basepars <- function(x, newdata, event){ 
  k <- event
  kpars <- (x$res$component == x$evnames[k])  & (x$res$dist == x$dists[k])
  bpars <- kpars & (x$res$terms %in% x$dlists[[k]]$pars)
  beta <- if (x$ncoveffsl[k]==0) 0 else x$res[kpars,"est"][x$covparsl[[k]]]
  if (is.null(newdata) | (x$ncoveffsl[k]==0)) {
    basepars <- as.list(x$res[bpars,"est"])
    names(basepars) <- x$dlists[[k]]$pars
    if (!is.null(newdata)) # if covariates for one event but not another, stretch out
      basepars <- lapply(basepars,function(x)x[rep(1,length(newdata[[1]])),drop=FALSE])
  }
  else { 
    basepars  <- vector(x$nthetal[k],  mode="list")
    names(basepars) <- x$dlists[[k]]$pars
    for (j in seq_len(x$nthetal[k])){
      parname <- x$dlists[[k]]$pars[j]
      betainds <- x$res$component == x$evnames[k] & x$res$parcov == x$dlists[[k]]$pars[j] 
      if (any(betainds)) { 
        tm <- delete.response(terms(x$all.formulae[[k]][[parname]]))
        if (!all(rownames(attr(tm, "factors")) %in% names(newdata))) 
          stop(sprintf("not all required variables supplied in `newdata`")) 
        xlev <- lapply(x$data$mf[[k]], levels)[attr(tm,"term.labels")]
        
        ## convert numeric newdata to factor if necessary, ensure ordered status preserved
        for (i in names(newdata)){
          vari <- x$data$mf[[1]][[i]]
          if (is.factor(vari)) {
            newdata[[i]] <- factor(as.character(newdata[[i]]), 
                                   levels=levels(vari), ordered=is.ordered(vari))
          }
        }

        mf <- model.frame(tm, newdata, xlev = xlev)
        mm <- model.matrix(tm, mf)
        basepar <- x$res[bpars,"est.t"][j]
        beta <-  x$res$est.t[betainds]
        Xb  <- as.numeric(mm %*% c(basepar, beta))
        basepars[[j]] <- x$dlists[[k]]$inv.transform[[j]](Xb)
      } else basepars[[j]] <- x$res[bpars,"est"][j]
    }
  }
  basepars
}

get_probpars <- function(x, newdata=NULL, na.action){ 
  if (is.null(newdata) | (x$ncoveffsp == 0) ) {
    probpars <- matrix(x$res$est[x$res$baseorcov=="pbase"], nrow=1)
  }
  else { 
    tm <- delete.response(terms(x$pformula))
    if (!all(rownames(attr(tm, "factors")) %in% names(newdata))) 
      stop(sprintf("not all required variables supplied in `newdata`")) 
    
    xlev <- lapply(x$data$mf[[1]], levels)[attr(tm,"term.labels")]
    mf <- model.frame(tm, newdata, xlev = xlev)
    X <- model.matrix(tm, mf)
    K <- x$K
    probpars <- matrix(nrow=nrow(X), ncol=K)
    for (i in 1:nrow(X)) {
      logitpcov <- numeric(K-1)
      for (k in 2:K){ 
        logitpbase <- x$res$est.t[grep(sprintf("prob%s$",k), x$res$terms)]
        beta <- x$res$est[grep(sprintf("prob%s\\(.+\\)",k), x$res$terms)]
        logitpcov[k-1] <-  X[i,,drop=FALSE] %*% c(logitpbase, beta)
      }
      probpars[i,] <- c(1, exp(logitpcov)) / (1  + sum(exp(logitpcov)))
    }
  }
  probpars <- as.data.frame(probpars)
  names(probpars) <- x$evnames
  if (!is.null(newdata))
    probpars <- cbind(probpars, newdata)
  probpars <- tidyr::pivot_longer(probpars, tidyselect::all_of(x$evnames), names_to="event", values_to="val")
  probpars 
}



##' Aalen-Johansen nonparametric estimates comparable to a fitted flexsurvmix
##' model
##'
##' Given a fitted flexsurvmix model, return the Aalen-Johansen estimates of the
##' probability of occupying each state at a series of times covering the
##' observed data.  State 1 represents not having experienced any of the
##' competing events, while state 2 and any further states correspond to having
##' experienced each of the competing events respectively.  These estimates can
##' be compared with the fitted probabilities returned by
##' \code{\link{p_flexsurvmix}} to check the fit of a \code{flexsurvmix} model.
##'
##' This is only supported for models with no covariates or models containing
##' only factor covariates.
##'
##' For models with factor covariates, the Aalen-Johansen estimates are computed
##' for the subsets of the data defined in \code{newdata}.  If \code{newdata} is
##' not supplied, then this function returns state occupancy probabilities for
##' all possible combinations of the factor levels.
##'
##' The Aalen-Johansen estimates are computed using
##' \code{\link[survival]{survfit}} from the \code{survival} package (Therneau
##' 2020).
##'
##' @param x Fitted model returned by \code{\link{flexsurvmix}}.
##'
##' @param newdata Data frame of alternative covariate values to check fit for.
##' Only factor covariates are supported. 
##'
##' @param tidy If \code{TRUE} then a single tidy data frame is returned.
##'   Otherwise the function returns the object returned by \code{survfit}, or a
##'   list of these objects if we are computing subset-specific estimates.
##'
##' @references Therneau T (2020). _A Package for Survival Analysis in R_. R
##'   package version 3.2-3, <URL: https://CRAN.R-project.org/package=survival>.
##'
##' @export
ajfit <- function(x, newdata=NULL, tidy=TRUE){ 
  dat <- x$data$mfcomb
  covnames <- attr(dat, "covnames")
  ncovs <- length(covnames)
  if (ncovs == 0){
    sf <- sftidy <- ajfit.dat(x$data$mf[[1]], x$evnames)
    ret <- if (tidy) as.data.frame(unclass(sf)[c("time","pstate","lower","upper")]) else sf
  }
  else { 
    faccovs <- sapply(dat[,covnames], is.factor)
    if (!all(faccovs)) 
      stop("Nonparametric estimation not supported with non-factor covariates")
    if (is.null(newdata)) {
      newdata <- do.call(expand.grid, lapply(dat[,covnames,drop=FALSE],  levels))
      for (i in names(newdata)) 
        class(newdata[[i]]) <- class(x$data$mf[[1]][[i]]) # preserve e.g. ordered factor status
    }
    newdata <- as.data.frame(newdata)
    ncovvals <- nrow(newdata)
    sf <- sftidy <- vector(ncovvals, mode="list")
    covnames <- names(newdata)
    ## TODO error handling 
    for (i in 1:ncovvals) { 
      datsub <- x$data$mf[[1]] # assume any of the mf[[k]] will do
      for (j in seq_along(covnames))  
        datsub <- datsub[datsub[,covnames[j]] == newdata[i,covnames[j]],,drop=FALSE]
      sf[[i]] <- ajfit.dat(datsub, x$evnames)
      sftidy[[i]] <- as.data.frame(unclass(sf[[i]])[c("time","pstate","lower","upper")])
      covvals <- newdata[i,,drop=FALSE][rep(1,nrow(sftidy[[i]])),,drop=FALSE]
      sftidy[[i]] <- cbind(sftidy[[i]], covvals)
    }
    ret <- if (tidy) do.call("rbind", sftidy) else sf
  }
  ret
}

ajfit.dat <- function(dat,evnames){
  intcens <-  model.response(dat)[,"status"] == 3
  dat <- dat[!intcens,]  ## exclude interval-censored data
  event <- model.extract(dat, "event")
  # create a multistate outcome for survfit.
  event <- as.character(event)
  event[is.na(event)] <- 0
  event <- factor(event, levels=c(0,evnames))
  Y <- model.response(dat)
  tname <- if ("time1" %in% colnames(Y)) "time1" else "time"
  time <- as.numeric(Y[,tname])
  sf <- survival::survfit(Surv(time, event) ~ 1)  
  sf
}



##' Forms a tidy data frame for plotting the fit of parametric mixture
##' multi-state models against nonparametric estimates
##'
##' This computes Aalen-Johansen estimates of the probability of occupying each
##' state at a series of times, using \code{\link{ajfit}}. The equivalent
##' estimates from the parametric model are then produced using
##' \code{\link{p_flexsurvmix}}, and concatenated with the nonparametric
##' estimates to form a tidy data frame. This data frame can then simply be
##' plotted using \code{\link[ggplot2]{ggplot}}.
##'
##' @param x Fitted model returned by \code{\link{flexsurvmix}}.
##'
##' @param maxt Maximum time to produce parametric estimates for.  By default
##'   this is the maximum event time in the data, the maximum time we have
##'   nonparametric estimates for.
##'
##' @param startname Label to give the state corresponding to "no event happened
##'   yet".  By default this is \code{"Start"}.
##'   
##' @param B Number of simulation replications to use to calculate a confidence
##' interval for the parametric estimates in \code{\link{p_flexsurvmix}}. 
##' Comparable intervals for the Aalen-Johansen estimates are returned if this
##' is set.  Otherwise if \code{B=NULL} then no intervals are returned.
##' 
##' @export
ajfit_flexsurvmix <- function(x, maxt=NULL, startname="Start", B=NULL){
  nstates <- x$K + 1
  dat <- x$data$mfcomb
  covnames <- attr(dat, "covnames")
  ncovs <- length(covnames)
  if (ncovs > 0){
    faccovs <- sapply(dat[,covnames,drop=FALSE], is.factor)
    if (!all(faccovs)) 
      stop("Nonparametric estimation not supported with non-factor covariates")
    newdata <- do.call(expand.grid, lapply(dat[,covnames,drop=FALSE],  levels))
    for (i in names(newdata)) 
      class(newdata[[i]]) <- class(dat[[i]]) # preserve e.g. ordered factor status
  } else newdata <- NULL    
  statenames <- c(startname,x$evnames)
  ajlong <- ajfit(x) %>% 
    tidyr::pivot_longer(cols = c(tidyselect::num_range("pstate.",1:nstates), 
                          tidyselect::num_range("lower.",1:nstates),
                          tidyselect::num_range("upper.",1:nstates)),
                 names_to=c("summary","state"),
                 names_sep="\\.", values_to="prob") %>%
    tidyr::pivot_wider(names_from="summary", values_from="prob")
  ajlong$state <- as.character(factor(ajlong$state, labels=statenames))
  names(ajlong)[names(ajlong)=="pstate"] <- "val"
  ajlong$model <- "Aalen-Johansen"
  if (is.null(B)) {
    ajlong$lower <- ajlong$upper <- NULL
    vals <- "val"
  } else {
    vals <- c("val","lower","upper")
  }
  names(ajlong)[names(ajlong)=="time"] <- "t"
  if (is.null(maxt)) maxt <- max(ajlong$t)
  times <- seq(0, maxt, length=100)
  ## FIXME p_flex returns factor, ajlong$x has ordered
  ## p_flex should return factor if that's what newdata has. so the newdata used should have ordered
  modcomp <- 
    p_flexsurvmix(x, t=times, newdata=newdata, startname=startname, B=B) %>%
    dplyr::mutate(model="Parametric mixture") %>%
    dplyr::full_join(ajlong, by = c("t", "state", names(newdata), "model", vals)) %>%
    dplyr::rename(time="t")
  modcomp
}
