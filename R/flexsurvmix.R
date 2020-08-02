##' Flexible parametric mixture models for times to competing events
##'
##' In a mixture model for competing events, an individual can experience one of
##' a set of different events.  We specify a model for the probability that they
##' will experience each event before the others, and a model for the time to
##' the event conditionally on that event occurring first.
##'
##' This differs from the more usual "competing risks" models, where we specify
##' "cause-specific hazards" describing the time to each competing event.  This
##' time will not be observed for an individual if one of the competing events
##' happens first.  The event that happens first is defined by the minimum of
##' the times to the alternative events.
##'
##' The \code{flexsurvmix} function fits a mixture model to data consisting of a
##' single time to an event for each individual, and an indicator for what type
##' of event occurs for that individual.   The time to event may be observed or
##' censored, just as in \code{\link{flexsurvreg}}, and the type of event may be
##' known or unknown. In a typical application, where we follow up a set of
##' individuals until they experience an event or a maximum follow-up time is
##' reached, the event type is known if the time is observed, and the event type
##' is unknown when follow-up ends and the time is right-censored.
##'
##' The model is fitted by maximum likelihood, either directly or by using an
##' expectation-maximisation (EM) algorithm, by wrapping
##' \code{\link{flexsurvreg}} to compute the likelihood or to implement the E
##' and M steps.
##'
##'
##' @param formula Survival model formula.  The left hand side is a \code{Surv}
##'   object specified as in \code{\link{flexsurvreg}}.  Any covariates on the
##'   right hand side of this formula will be placed on the location parameter
##'   for every component-specific distribution.  Different covariates may be
##'   supplied on different components using the \code{anc} argument.
##'
##' @param data Data frame containing variables mentioned in \code{formula},
##'   \code{event} and \code{anc}.
##'
##' @param event Variable in the data that specifies which of the alternative
##'   events is observed for which individual.  If the individual's follow-up is
##'   right-censored, or if the event is otherwise unknown, this variable must
##'   have the value \code{NA}.
##'
##'   Ideally this should be a factor, since the mixture components can then be
##'   easily identified in the results with a name instead of a number.  If this
##'   is not already a factor, it is coerced to one.   Then the levels of the
##'   factor define the required order for the components of the list arguments
##'   \code{dists}, \code{anc}, \code{inits} and \code{dfns}.  Alternatively, if
##'   the components of the list arguments are named according to the levels of
##'   \code{event}, then the components can be arranged in any order.
##'
##' @param dists Vector specifying the parametric distribution to use for each
##'   component. The same distributions are supported as in
##'   \code{\link{flexsurvreg}}.
##'
##' @param pformula Formula describing covariates to include on the component
##'   membership proabilities by multinomial logistic regression.  The first
##'   component is treated as the baseline.
##'
##' @param anc List of component-specific lists, of length equal to the number
##'   of components.   Each component-specific list is a list of formulae
##'   representing covariate effects on parameters of the distribution.
##'
##'   If there are covariates for one component but not others, then a list
##'   containing one null formula on the location parameter should be supplied
##'   for the component with no covariates, e.g \code{list(rate=~1)} if the
##'   location parameter is called \code{rate}.
##'
##'   Covariates on the location parameter may also be supplied here instead of
##'   in \code{formula}.  Supplying them in \code{anc} allows some components
##'   but not others to have covariates on their location parameter.
##'
##' @param initp Initial values for component membership probabilities.  By
##'   default, these are assumed to be equal for each component.
##'
##' @param inits List of component-specific vectors. Each component-specific
##'   vector contains the initial values for the parameters of the
##'   component-specific model, as would be supplied to
##'   \code{\link{flexsurvreg}}.   By default, a heuristic is used to obtain
##'   initial values, which depends on the parametric distribution being used,
##'   but is usually based on the empirical mean and/or variance of the survial
##'   times.
##'
##' @param fixedpars Indexes of parameters to fix at their initial values and
##'   not optimise. Arranged in the order: baseline mixing probabilities,
##'   covariates on mixing probabilities, time-to-event parameters by mixing
##'   component.  Within mixing components, time-to-event parameters are ordered
##'   in the same way as in \code{\link{flexsurvreg}}.
##'
##'   If \code{fixedpars=TRUE} then all parameters will be fixed and the
##'   function simply calculates the log-likelihood at the initial values.
##'
##'   Not currently supported when using the EM algorithm.
##'
##' @param dfns List of lists of user-defined distribution functions, one for
##'   each mixture component.  Each list component is specified as the
##'   \code{dfns} argument of \code{\link{flexsurvreg}}.
##'
##' @param method Method for maximising the likelihood.  Either \code{"em"} for
##'   the EM algorithm, or \code{"direct"} for direct maximisation.  EM is only
##'   currenly supported if there are no covariates.
##'
##' @param em.control List of settings to control EM algorithm fitting.  The
##'   only options currently available are
##'
##'   \code{trace} set to 1 to print the parameter estimates at each iteration
##'   of the EM algorithm
##'
##'   \code{reltol} convergence criterion.  The algorithm stops if the log
##'   likelihood changes by a relative amount less than \code{reltol}.  The
##'   default is the same as in \code{\link{optim}}, that is,
##'   \code{sqrt(.Machine$double.eps)}.
##'
##'   For example, \code{em.control = list(trace=1, reltol=1e-12)}.
##'
##' @param optim.control List of options to pass as the \code{control} argument
##'   to \code{\link{optim}},  which is used by \code{method="direct"} or in the
##'   M step of \code{method="em"}.  By default, this uses \code{fnscale=10000}
##'   and \code{ndeps=rep(1e-06,p)} where \code{p} is the number of parameters
##'   being estimated, unless the user specifies these options explicitly.
##'
##'   

##' @inheritParams flexsurvreg
##'
##' @return List of objects containing information about the fitted model.   The
##'   important one is \code{res}, a data frame containing the parameter
##'   estimates and associated information.
##'
##'
##' @references Larson, M. G., & Dinse, G. E. (1985). A mixture model for the
##' regression analysis of competing risks data. Journal of the Royal
##' Statistical Society: Series C (Applied Statistics), 34(3), 201-211.
##'
##' Lau, B., Cole, S. R., & Gange, S. J. (2009). Competing risk regression
##' models for epidemiologic data. American Journal of Epidemiology, 170(2),
##' 244-256.
##'
##' @export
flexsurvmix <- function(formula, data, event, dists,
                        pformula=NULL, anc=NULL, 
                        initp=NULL, inits=NULL, 
                        fixedpars=NULL, dfns=NULL,
                        method="direct", 
                        em.control=NULL, 
                        optim.control=NULL,
                        aux=NULL, 
                        sr.control=survreg.control(), 
                        integ.opts, ...){
  call <- match.call()
  ## Determine names of competing events, and their order, based on event data 
  event <- eval(substitute(event), data, parent.frame())
  if (!is.factor(event)) event <- factor(event)
  evnames <- levels(event)
  
  ## For all arguments that are vectors or lists by event, 
  ## make sure their names match names of events
  ## and reorder if necessary so the order of the names matches too
  if (missing(dists)) stop("Distributions \"dists\" not specified")
  dists <- clean_listarg(dists, "dists", evnames)
  anc <- clean_listarg(anc, "anc", evnames)
  inits <- clean_listarg(inits, "inits", evnames)
  dfns <- clean_listarg(dfns, "dfns", evnames)
  
  ## Build distribution functions
  K <- length(dists)
  dlists <- vector(K, mode="list")
  if (is.null(dfns)) dfns <- vector(K, mode="list")
  for (k in 1:K){
    dlists[[k]] <- parse.dist(dists[k])
    dfns[[k]] <- form.dp(dlists[[k]], dfns[[k]], integ.opts)
  }
  names(dlists) <- names(dfns) <- names(dists)
  check.formula(formula, dlists[[1]])
  
  ## Build covariate model formulae
  if (!is.null(anc)) {
    if (!is.list(anc)) stop("`anc` should be a list")
  }
  ancm <- vector(K, mode="list")
  for (k in 1:K){
    msg <- sprintf("anc[[%s]] must be a list of formulae", k)
    ancm[[k]] <- anc_from_formula(formula, anc[[k]], dlists[[k]], msg)
  }
  locform <- forms <- vector(K, mode="list") 
  for (k in 1:K) { 
    ancnames <- setdiff(dlists[[k]]$pars, dlists[[k]]$location)
    locform[[k]]  <- get.locform(formula, ancnames)
    loc <- dlists[[k]]$location
    if (loc %in% names(ancm[[k]])){
      locform[[k]] <- update(locform[[k]], ancm[[k]][[loc]])
      ancm[[k]][[loc]] <- anc[[k]][[loc]]<- NULL
    }
    forms[[k]] <- c(location=locform[[k]], ancm[[k]])
    names(forms[[k]])[[1]] <- loc
  }
  f2 <- concat.formulae(formula, c(unlist(forms),pformula))
  
  ## Build model frame given formulae
  indx <- match(c("formula", "data", "event"), names(call), nomatch = 0)
  if (indx[1] == 0)
    stop("A \"formula\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  temp[["event"]] <- event
  temp[["formula"]] <- f2
  if (missing(data)) temp[["data"]] <- environment(formula)
  temp[["na.action"]] <- na.pass # event will have NAs by design. Or recode them ? TESTME
  m <- eval(temp, parent.frame())
  m <- droplevels(m) # remove unused factor levels after subset applied
  nobs <- nrow(m) 
  
  ## Build design matrices given formulae and model frame  
  mml <- mx <- X <- whichparcov <- vector(K, mode="list")
  for (k in 1:K) { 
    mml[[k]] <- mx[[k]] <- vector(mode="list", length=length(dlists[[k]]$pars))
    names(mml[[k]]) <- names(mx[[k]]) <- c(dlists[[k]]$location, 
                                           setdiff(dlists[[k]]$pars, dlists[[k]]$location))
    for (i in names(forms[[k]])){
      mml[[k]][[i]] <- model.matrix(forms[[k]][[i]], m)
      mx[[k]][[i]] <- length(unlist(mx[[k]])) + seq_len(ncol(mml[[k]][[i]][,-1,drop=FALSE]))
    }
    X[[k]] <- compress.model.matrices(mml[[k]])
    npc <- unlist(lapply(mml[[k]], ncol))  - 1
    whichparcov[[k]] <- rep(names(npc), npc)
  }
  if (is.null(pformula)) pformula <- ~1
  Xp <- model.matrix(pformula, m)[,-1,drop=FALSE]
  
  ## Convert event internally to numbers based on factor levels
  event <- model.extract(m, "event")
  event <- match(event, evnames, incomparables=NA)
  
  ## Count parameters of particular kinds 
  ncoveffsl <- sapply(X, ncol) # number of covariate effects for each component-specific distribution
  ncoveffs <- sum(ncoveffsl)    
  ncovsp <- ncol(Xp)  # number of covariates on mixing probabilities
  ncoveffsp <- (K-1)*ncovsp  # number of covariate effects on mixing probabilities
  nppars <- length(dists) - 1 + ncoveffsp # number of parameters related to mixing probabilities: baseline probs and covariate effects
  ## Order of time-to-event parameters: baseline parameters in the order they
  ## appear in the distribution function  (e.g. shape before scale in the
  ## Weibull), followed by covariate effects on the location parameter, followed
  ## by covariate effects on the remaining parameters.
  nthetal <- sapply(dlists, function(x)length(x$pars)) # number of baseline parameters for each component-specific distribution
  ncparsl <- nthetal + ncoveffsl # total number of pars for each component
  parindsl <- rep(1:K, ncparsl) # index identifying component for each of those parameters 
  ncpars <- sum(ncparsl) # total number of parameters related to time-to-event distributions 
  npars <- nppars + ncpars # total number of parameters 
  ## identify parameters of particular kinds 
  parclass <- rep(c("prob","time"), c(nppars, ncpars))  # mixing prob or time-to-event related 
  baseorcov <- c(rep(c("pbase","pcov"), c(length(dists) -1, ncoveffsp)),  # baseline or covariate effect 
                 rep(rep(c("tbase","tcov"), K), as.vector(rbind(nthetal, ncoveffsl))))
  parcov <- character(npars) # identify baseline parameter which a covariate effect modifies
  parcov[parclass=="prob"  &  baseorcov=="pcov"] <- rep(paste0("prob",2:K), each=ncovsp)
  parcov[baseorcov=="tcov"] <- unlist(whichparcov)
  
  ## Build initial values for each parameter type for optimisation 
  if (is.null(initp)) initp <- rep(1/K, K)
  alpha <- initp 
  Y <- check.flexsurv.response(model.extract(m, "response"))
  theta_inits <- cov_inits <- covparsl <- vector(K, mode="list")
  for (k in 1:K) { 
    covparsl[[k]] <- nthetal[k] + seq_len(ncoveffsl[k])
    if (is.null(inits[[k]])) { 
      yy <- ifelse(Y[,"status"]==3 & is.finite(Y[,"time2"]), (Y[,"time1"] + Y[,"time2"])/2, Y[,"time1"])
      yy  <- yy[event==k | is.na(event)]
      # wt <- yy*weights*length(yy)/sum(weights)
      dlists[[k]]$inits <- expand.inits.args(dlists[[k]]$inits)
      ## This extra stuff is needed for the Weibull, to use a  survreg fit to
      ## get "initial values"
      inits.aux <- c(aux, list(forms=forms[[k]], data=if(missing(data)) NULL else data, weights=temp$weights,
                               control=sr.control,
                               counting=(attr(model.extract(m, "response"), "type")=="counting")
      ))
      ## Auto-generate initial values using the heuristic for that distribution
      auto.inits <- dlists[[k]]$inits(t=yy,mf=m,mml=mml[[k]],aux=inits.aux)
      nin <- length(inits)
      theta_inits[[k]] <- auto.inits[seq_len(nthetal[k])]
      if ((length(auto.inits) > nthetal[k]) && (ncoveffsl[k] > 0)) { 
        ## e.g for Weibull distributions, initial value fn also estimates covariate effects. 
        cov_inits[[k]] <- auto.inits[nthetal[k] + seq_len(ncoveffsl[k])]
      } else 
        cov_inits[[k]] <- rep(0, ncoveffsl[k])
    } else { 
      theta_inits[[k]] <- inits[[k]][1:nthetal[k]] # baseline pars of parametric survival dists
      cov_inits[[k]] <- inits[[k]][covparsl[[k]]]
    }
    names(theta_inits[[k]]) <- dlists[[k]]$pars
    names(cov_inits[[k]]) <- colnames(X[[k]])
  }
  ## Transform initial values to log or logit scale for optimisation, if needed. 
  inits_probs <- initp
  inits_alpha <- qlogis(inits_probs[2:K])
  inits_covp <- rep( rep(0,ncovsp), K-1)   # order by covariate within probability
  names(inits_covp) <- sprintf("prob%s(%s)", rep(2:K, each=ncovsp), rep(colnames(Xp), K-1))
  inits_theta <- numeric()
  for (k in 1:K){
    inits_theta <- c(inits_theta, par.transform(theta_inits[[k]], dlists[[k]]), cov_inits[[k]])
  }
  
  loglik_flexsurvmix <- function(parsopt, ...){
    pars <- inits_all
    pars[optpars] <- parsopt
    ## Apply covariate effects to mixing probabilities by multinomial logistic
    ## regression
    alpha <-  c(0, pars[1:(K-1)])
    alphamat <- matrix(alpha, nrow=nobs, ncol=K, byrow=TRUE)  # by individual 
    if (ncovsp > 0) { 
      for (k in 2:K){  
        cpinds <- K - 1 + (k-2)*ncovsp + 1:ncovsp
        alphamat[,k] <- alpha[k] + Xp %*% pars[cpinds] 
      }
    }
    pmat <- exp(alphamat)
    pmat <- pmat / rowSums(pmat)
    probmat <- pmat # this will be 1 or 0 if event observed
    
    ## Remaining parameters are time-to-event stuff   
    parsl <- split(pars[nppars + seq_len(ncpars)], parindsl)
    theta <- coveffs <- vector(K, mode="list")
    
    ## Contribution to likelihood for mixing probabilities for those with known events
    llp_event_known <- matrix(0, nrow=nobs, ncol=K)
    
    ## Likelihood from times to events or censoring
    liki <-  matrix(0, nrow=nobs, ncol=K)
    for (k in 1:K){
      ## Evaluate likelihood for each component by calling flexsurvreg with
      ## parameters fixed at initial values. Build this initial values vector.
      theta[[k]] <- inv.transform(parsl[[k]][seq_len(nthetal[k])], dlists[[k]])
      coveffs[[k]] <- parsl[[k]][nthetal[k] + seq_len(ncoveffsl[k])]
      initsk <- c(theta[[k]], coveffs[[k]])
      if (any(is.infinite(log(exp(initsk))))) return(-Inf)
      ## Event probability is 1 if their event was observed to be k
      probmat[!is.na(event) & event==k, k] <- 1  
      ## Event probability is 1 if their event was observed to be one other than k
      probmat[!is.na(event) & event!=k, k] <- 0
      need_lik <- probmat[,k] > 0
      liki[need_lik,k] <- exp(do.call("flexsurvreg", 
                                      list(formula=locform[[k]], data=data, dist=dists[k], 
                                           anc=anc[[k]], inits=initsk, subset=need_lik,
                                           fixedpars=TRUE))$logliki)
      llp_event_known[,k] <- as.numeric((!is.na(event) & event==k) * log(pmat[,k]))
    }
    logliki <- - (log(rowSums(probmat*liki, na.rm=TRUE)) + rowSums(llp_event_known))
    res <- sum(logliki) 
    attr(res, "indiv") <- logliki
    res
  }
  # Parameter dictionary to be completed with the estimates
  res <- data.frame(component = c(evnames, rep(evnames[1:(K-1)], each=ncovsp),  rep(evnames, ncparsl)), 
                    dist = c(rep("",nppars+1), rep(dists, ncparsl)), 
                    terms = c(paste0("prob",1:K),  names(inits_covp),  names(inits_theta)),
                    parclass = c("prob", parclass), 
                    baseorcov = c("pbase", baseorcov), 
                    parcov = c("", parcov))
  inits_all <- c(inits_alpha, inits_covp, inits_theta)
  if (isTRUE(fixedpars)) fixedpars <- seq_len(npars)
  optpars <- setdiff(seq_len(npars), fixedpars)
  fixed <- !any(optpars)
  inits_opt <- inits_all[optpars]
  if (any(fixedpars) && (method=="em")){
    method <- "direct"
    if (any(optpars)) 
      warning("Optimisation with some parameters fixed not currently supported with EM algorithm, switching to direct likelihood maximisation")
  }
  if (is.null(optim.control$fnscale))
    optim.control$fnscale <- 10000
  
  if (method=="direct"){
    if (any(optpars)){
      if (is.null(optim.control$ndeps)) 
        optim.control$ndeps = rep(1e-06, length(inits_opt))
      opt <- optim(inits_opt, loglik_flexsurvmix, hessian=TRUE, method="BFGS",
                   control=optim.control, ...)
    } else {
      opt <- list(par=inits_all, value=loglik_flexsurvmix(inits_all))
    }
    logliki <- -attr(loglik_flexsurvmix(opt$par), "indiv")
    
    ## Transform mixing probs back to natural scale for presentation 
    opt_all <- inits_all
    opt_all[optpars] <- opt$par
    res_alpha <- plogis(opt_all[1:(K-1)])
    res_alpha <- c(1 - sum(res_alpha), res_alpha)
    res_covp <- if (ncoveffsp>0) opt_all[K:(K-1+ncoveffsp)] else numeric()
    res_cpars <- split(opt_all[(nppars+1):length(opt_all)], parindsl)
    res_theta <- res_coveffs <- vector(K, mode="list")
    for (k in 1:K){
      ## Transform time-to-event parameters back to natural scale 
      res_theta[[k]] <- inv.transform(res_cpars[[k]][seq_len(nthetal[k])], dlists[[k]])
      res_coveffs[[k]] <- res_cpars[[k]][nthetal[k] + seq_len(ncoveffsl[k])]
      res_cpars[[k]] <- c(res_theta[[k]], res_coveffs[[k]])
    }
    ## Build tidy data frame of results with one row per parameter. Includes an
    ## extra row for the probability of the first event (defined as 1 - sum of
    ## rest)
    res <- cbind(res, data.frame(
                      est = c(res_alpha, res_covp, unlist(res_cpars)),
                      est.t = c(NA, opt_all)))
                      
    loglik <- - opt$value
    if (!fixed)
      cov <- solve(opt$hessian)
  }
  
  else if (method=="em") { 
    if (!is.list(em.control)) em.control <- list()
    if (is.null(em.control$reltol)) em.control$reltol <- sqrt(.Machine$double.eps)
    converged <- FALSE
    iter <- 0
    alpha <- inits_alpha
    covp <- inits_covp 
    theta <- theta_inits
    
    while (!converged) {  
      ## E step
      alphamat <- matrix(c(0,alpha), nrow=nobs, ncol=K, byrow=TRUE)  # by individual 
      if (ncovsp>0){
        for (k in 2:K){ 
          cpinds <- (k-2)*ncovsp + 1:ncovsp
          alphamat[,k] <- alpha[k-1] + Xp %*% covp[cpinds] 
        }
      }
      pmat <- exp(alphamat)
      pmat <- pmat / rowSums(pmat)
      
      llmat <- matrix(nrow=nobs, ncol=K)
      for (k in 1:K) {
        fs <- flexsurvreg(formula=formula, data=data, dist=dists[k], anc=anc[[k]], inits=theta[[k]], fixedpars=TRUE) 
        llmat[,k] <-  fs$logliki
      }
      alphap <- exp(llmat) * pmat
      w <- alphap  / rowSums(alphap) # probability that each observation belongs to each component
      for (k in 1:K) {
        w[!is.na(event) & event==k, k] <- 1
        w[!is.na(event) & event!=k, k] <- 0
      }
      w
      
      ## M step
      
      ## estimate component probs and covariate effects on them
      #alpha <- colMeans(w)
      loglik_p_em <- function(pars){
        alphamat <- matrix(c(0,pars[1:(K-1)]), nrow=nobs, ncol=K, byrow=TRUE)
        if (ncovsp > 0) { 
          for (k in 2:K){  
            cpinds <- K - 1 + (k-2)*ncovsp + 1:ncovsp
            alphamat[,k] <- alphamat[,k-1] + Xp %*% pars[cpinds] 
          }
        }
        pmat <- exp(alphamat)
        pmat <- pmat / rowSums(pmat)
        -sum(w*log(pmat))
      }

      #  if setting a control argument here in  the future, note ndeps has  to be different  
      # length from others. 
      parsp <- optim(c(alpha,covp), loglik_p_em, method="BFGS")
      alpha <- parsp$par[1:(K-1)]
      probs <- c(1-sum(plogis(alpha)), plogis(alpha))
      covp <- if (ncovsp>0) parsp$par[K:nppars] else NULL
      
      ## call flexsurvreg for each component on weighted dataset to estimate component-specific pars 
      thetanew <- vector(K, mode="list")
      ll <- numeric(K)
      
      for (k in 1:K) {
        if (is.null(optim.control$ndeps)) 
          optim.control$ndeps = rep(1e-06, length(theta[[k]]))
        fs <- do.call("flexsurvreg",  # need do.call to avoid environment faff with supplying weights
                      list(formula=formula, data=data, dist=dists[k], 
                           anc=anc[[k]], inits=theta[[k]], weights=w[,k], 
                           subset=w[,k]>0, hessian=FALSE,
                           control=optim.control))
        ll[k] <- fs$loglik
        thetanew[[k]] <- fs$res[,"est"]
      }
      
      theta <- thetanew
      logliknew <- sum(ll)
      if (iter > 0)
        converged <-  (abs(logliknew / loglik - 1) <= em.control$reltol)
      loglik <- logliknew
      est <- c(probs, covp, unlist(theta))
      est.t <- c(NA, alpha, covp, unlist(parlist.transform(theta,dlists)))
      if (is.numeric(em.control$trace) && em.control$trace > 0)
        print(est)
      iter <- iter + 1
    }
    res <- cbind(res, 
                 data.frame(est=est, est.t=est.t))
    ll <- loglik_flexsurvmix(est.t[-1])
    loglik <- -as.numeric(ll)
    logliki <- -attr(ll, "indiv")
    cov <- solve(numDeriv::hessian(loglik_flexsurvmix, est.t[-1]))
    opt <- NULL
  }
  
  ## Add standard errors to results data frame, given covariance matrix.
  ## var ( 1 - p1 - p2 - .. ) = var(p1) + var(p2) - cov(p1,p2) - ...
  res$fixed <- c(NA, rep(FALSE, npars))
  res$fixed[1 + fixedpars] <- TRUE
  if (!fixed){
    covp <- cov[1:(K-1), 1:(K-1), drop=FALSE]
    sepk <- sqrt(sum(diag(covp)) - sum(covp[lower.tri(covp)]))
    res$se <- rep(NA, npars+1)
    res$se[1] <- sepk
    res$se[1 + optpars] <- sqrt(diag(cov))
  }

  res <- list(call=match.call(),
              res=res, loglik=loglik,  cov=cov, 
              npars=npars, AIC=-2*loglik + 2*npars,
              K=K, dists=dists,  dlists=dlists, dfns=dfns, 
              opt=opt,
              fixedpars=fixedpars, optpars=optpars,
              logliki=logliki,
              evnames=evnames,
              nthetal=nthetal, parindsl=parindsl, 
              ncoveffsl=ncoveffsl,ncoveffsp=ncoveffsp, covparsl=covparsl,
              all.formulae=forms, pformula=pformula, 
              data=list(mf=m), mx=mx)
  class(res) <- "flexsurvmix"
  res
}

logLik.flexsurvmix <- function(object, ...){
  val <- object$loglik
  attributes(val) <- NULL
  attr(val, "df") <- object$npars
  attr(val, "nobs") <- nrow(model.frame(object))
  class(val) <- "logLik"
  val
}

clean_listarg <- function(arg, argname, evnames){
  if (!is.null(arg)){
    if (length(arg) != length(evnames))
      stop(sprintf("`%s` of length %s, should equal the number of mixture components, assumed to be %s, from the number of levels of factor(event)", 
           argname, length(arg), length(evnames)))
    narg <- names(arg)
    if (is.null(narg)) 
      names(arg) <- evnames
    else { 
      if (!identical(sort(narg), sort(evnames)))
        sprintf("names(%s) = %s, but this should match levels(factor(event)) = %s", 
                narg, argname, dput(evnames))
      arg <- arg[evnames]
    }
  }
  arg
}


##' @export
print.flexsurvmix <- function(x, ...)
{
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  if (x$npars > 0) {
    res <- x$res
    cat ("Estimates: \n")
    args <- list(...)
    if (is.null(args$digits)) args$digits <- 3
    f <- do.call("format", c(list(x=res), args))
    print(f, print.gap=2, quote=FALSE, na.print="")
  }
  cat("\nLog-likelihood = ", x$loglik, ", df = ", x$npars,
      "\nAIC = ", x$AIC, "\n\n", sep="")
}


par.transform <- function(pars, dlist){
  pars.t <- numeric(length(pars))
  names(pars.t) <- names(pars)
  for (i in seq_along(pars.t))
    pars.t[i] <- dlist$transforms[[i]](pars[i])
  pars.t
}

inv.transform <- function(pars, dlist){
  pars.nat <- numeric(length(pars))
  names(pars.nat) <- names(pars)
  for (i in seq_along(pars.nat))
    pars.nat[i] <- dlist$inv.transforms[[i]](pars[i])
  pars.nat
}

parlist.transform <- function(parlist, dlists){
  parlist.t <- vector(length(parlist), mode="list")
  for (k in seq_along(parlist.t))
    parlist.t[[k]] <- par.transform(parlist[[k]], dlists[[k]])
  parlist.t
}

invlist.transform <- function(parlist, dlists){
  parlist.t <- vector(length(parlist), mode="list")
  for (k in seq_along(parlist.t))
    parlist.t[[k]] <- inv.transform(parlist[[k]], dlists[[k]])  
  parlist.t
}

## Transform full vector of estimates back to natural scale.  
## Returns Kth baseline prob as well 
## TODO can this go into the main flexsurvreg function? Would have to form list first 

inv.transform.res <- function(x, dlists) {
  dpars <- x$res$est.t[x$res$dist!=""]
  dpars <- split(dpars, x$parindsl)
  dnames <- split(x$res$terms[x$res$dist!=""], x$parindsl)
  K <- x$K
  est <- vector(K, mode="list")
  for (k in 1:K) {
    bpars <- dpars[[k]][dnames[[k]] %in% x$dlists[[k]]$pars]
    cpars <- dpars[[k]][!(dnames[[k]] %in% x$dlists[[k]]$pars)]
    bpars <- inv.transform(bpars, dlists[[k]])
    est[[k]] <- c(bpars, cpars)
  }
  probs <- plogis(x$res$est.t[1:K])  # prob of group k given cov value of zero 
  probs[K] <- 1 - sum(probs[1:(K-1)]) 
  pcov <- x$res$est.t[grep("prob[[:digit:]]+\\(.+\\)", x$res$terms)]
  c(probs, pcov, unlist(est))
}

##' @export
model.frame.flexsurvmix <- function(formula, ...)
{
  x <- formula
  x$data$m
}
