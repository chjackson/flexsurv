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
##'   object specified as in \code{\link{flexsurvreg}}.  This may define various
##'   kinds of censoring, as described in \code{\link{Surv}}. Any covariates on
##'   the right hand side of this formula will be placed on the location
##'   parameter for every component-specific distribution. Covariates on other
##'   parameters of the component-specific distributions may be supplied  using
##'   the \code{anc} argument.
##'
##'   Alternatively, \code{formula} may be a list of formulae, with one
##'   component for each alternative event.  This may be used to specify
##'   different covariates on the location parameter for different components.
##'
##'   A list of formulae may also be used to indicate that for particular
##'   individuals, different events may be observed in different ways, with
##'   different censoring mechanisms.  Each  list component specifies the data
##'   and censoring scheme for that mixture component.
##'
##'   For example, suppose we are studying people admitted to hospital,and the
##'   competing states are death in hospital and discharge from hospital.  At
##'   time t we know that a particular individual is still alive, but we do not
##'   know whether they are still in hospital, or have been discharged.  In this
##'   case, if the individual were to die in hospital, their death time would be
##'   right censored at t.  If the individual will be (or has been) discharged
##'   before death, their discharge time is completely unknown, thus
##'   interval-censored on (0,Inf). Therefore,  we need to store different event
##'   time and status variables in the data for different alternative events.
##'   This is specified here as
##'
##'   \code{formula = list("discharge" = Surv(t1di, t2di, type="interval2"),
##'   "death" = Surv(t1de, status_de))}
##'
##'   where for this individual, \code{(t1di, t2di) = (0, Inf)} and \code{(t1de,
##'   status_de)  = (t, 0)}.
##'
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
##'   but not others to have covariates on their location parameter.  If a covariate
##'   on the location parameter was provided in \code{formula}, and there are 
##'   covariates on other parameters, then a null formula should be included 
##'   for the location parameter in \code{anc}, e.g \code{list(rate=~1)}
##'
##' @param partial_events List specifying the factor levels of \code{event}
##'   which indicate knowledge that an individual will not experience particular
##'   events, but may experience others.   The names of the list indicate codes
##'   that indicate partial knowledge for some individuals.  The list component
##'   is a vector, which must be a subset of \code{levels(event)} defining the
##'   events that a person with the corresponding event code may experience.
##'
##'   For example, suppose there are three alternative events called
##'   \code{"disease1"},\code{"disease2"} and \code{"disease3"}, and for some
##'   individuals we know that they will not experience \code{"disease2"}, but
##'   they may experience the other two events.  In that case we must create a
##'   new factor level, called, for example \code{"disease1or3"}, and set the
##'   value of \code{event} to be \code{"disease1or3"} for those individuals.
##'   Then we use the \code{"partial_events"} argument to tell
##'   \code{flexsurvmix} what the potential events are for individuals with this
##'   new factor level.
##'
##'   \code{partial_events = list("disease1or3" = c("disease1","disease3"))}
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
##'   the EM algorithm, or \code{"direct"} for direct maximisation.
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
##'   \code{var.method} method to compute the covariance matrix. \code{"louis"}
##'   for the method of Louis (1982), or \code{"direct"}for direct numerical
##'   calculation of the Hessian of the log likelihood.
##'
##'   \code{optim.p.control} A list that is passed as the \code{control}
##'   argument to  \code{optim} in the M step for the component membership
##'   probability parameters. The optimisation in the M step for the
##'   time-to-event parameters can be controlled by the \code{optim.control}
##'   argument to \code{flexsurvmix}.
##'
##'   For example, \code{em.control = list(trace=1, reltol=1e-12)}.
##'
##' @param optim.control List of options to pass as the \code{control} argument
##'   to \code{\link{optim}},  which is used by \code{method="direct"} or in the
##'   M step for the time-to-event parameters in \code{method="em"}.  By
##'   default, this uses \code{fnscale=10000} and \code{ndeps=rep(1e-06,p)}
##'   where \code{p} is the number of parameters being estimated, unless the
##'   user specifies these options explicitly.
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
##'   regression analysis of competing risks data. Journal of the Royal
##'   Statistical Society: Series C (Applied Statistics), 34(3), 201-211.
##'
##'   Lau, B., Cole, S. R., & Gange, S. J. (2009). Competing risk regression
##'   models for epidemiologic data. American Journal of Epidemiology, 170(2),
##'   244-256.
##'
##' @export
flexsurvmix <- function(formula, data, event, dists,
                        pformula=NULL, anc=NULL,
                        partial_events = NULL,
                        initp=NULL, inits=NULL,
                        fixedpars=NULL, 
                        dfns=NULL,
                        method="direct",
                        em.control=NULL,
                        optim.control=NULL,
                        aux=NULL,
                        sr.control=survreg.control(),
                        integ.opts, hess.control=NULL, ...){
  call <- match.call()
  ## Determine names of competing events, and their order, based on event data
  event <- eval(substitute(event), data, parent.frame())
  if (!is.factor(event)) event <- factor(event)

  ## Determine which of these are codes for partially observed events
  if (!is.null(partial_events)) {
    if (!is.list(partial_events)) stop("`partial_events` must be a list")
    if (!is.null(names(partial_events))) stop("`partial_events` must be a named list")
    npartials <- length(partial_events)
    evnames <- setdiff(levels(event), names(partial_events))
    for (i in 1:npartials){
      badnames <- partial_events[[i]][!(partial_events[[i]] %in% evnames)]
      badnames <- paste0("\"",badnames,"\"")
      if (length(badnames) > 0)
        stop(sprintf("partial_events[[%s]] values %s not in levels(events)",
                     i, paste(badnames,collapse=",")))
      ## Convert to codes
      partial_events[[i]] <- match(partial_events[[i]], evnames)
    }
  } else npartials <- 0
  evnames <- setdiff(levels(event), names(partial_events))

  partial_event <- ifelse(event %in% names(partial_events), event, NA)
  partial_event <- match(partial_event, names(partial_events))
  event[event %in% names(partial_events)] <- NA

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
    dlists[[k]] <- parse.dist(dists[[k]])
    dfns[[k]] <- form.dp(dlists[[k]], dfns[[k]], integ.opts)
  }
  names(dlists) <- names(dfns) <- names(dists)

  check.formula.flexsurvmix(formula, dlists[[1]]) # TESTME
  if (is.list(formula)) formula <- clean_listarg(formula, "formula", evnames)

  ## Build covariate model formulae
  if (!is.null(anc)) {
    if (!is.list(anc)) stop("`anc` should be a list")
  }
  ancm <- vector(K, mode="list")
  for (k in 1:K){
    msg <- sprintf("anc[[%s]] must be a list of formulae", k)
    fk <- if (is.list(formula)) formula[[k]] else formula
    ancm[[k]] <- anc_from_formula(fk, anc[[k]], dlists[[k]], msg)
  }
  locform <- forms <- vector(K, mode="list")
  for (k in 1:K) {
    ancnames <- setdiff(dlists[[k]]$pars, dlists[[k]]$location)
    fk <- if (is.list(formula)) formula[[k]] else formula
    locform[[k]]  <- get.locform(fk, ancnames)
    loc <- dlists[[k]]$location
    if (loc %in% names(ancm[[k]])){
      ## Add any extra covariates on the location parameter found in "anc" (usually we'll be adding to ~1)
      fc <- as.character(ancm[[k]][[loc]])
      extraterms <- formula(paste(fc[1],  ". +", fc[-1], collapse=" "))
      locform[[k]] <- update(locform[[k]], extraterms)
      ancm[[k]][[loc]] <- anc[[k]][[loc]]<- NULL
    }
    forms[[k]] <- c(location=locform[[k]], ancm[[k]])
    names(forms[[k]])[[1]] <- loc
  }

  ## Build model frame given formulaec
  indx <- match(c("formula", "data", "event"), names(call), nomatch = 0)
  if (indx[1] == 0)
    stop("A \"formula\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  temp[["event"]] <- event
  f1 <- if (is.list(formula)) formula[[1]] else formula
  if (missing(data)) temp[["data"]] <- environment(f1)
  temp[["na.action"]] <- na.pass # event will have NAs by design. Or recode them ? TESTME

  ## one model frame for each component.  may be different if different response terms in formula
  m <- vector(K, mode="list")
  for (k in 1:K){
    f2 <- concat.formulae(locform[[k]], c(unlist(forms), pformula))
    temp[["formula"]] <- f2
    m[[k]] <- eval(temp, parent.frame())
    m[[k]] <- droplevels(m[[k]]) # remove unused factor levels after subset applied
  }
  nobs <- nrow(m[[1]])

  ## Build design matrices given formulae and model frame
  mml <- mx <- X <- whichparcov <- vector(K, mode="list")
  for (k in 1:K) {
    mml[[k]] <- mx[[k]] <- vector(mode="list", length=length(dlists[[k]]$pars))
    names(mml[[k]]) <- names(mx[[k]]) <- c(dlists[[k]]$location,
                                           setdiff(dlists[[k]]$pars, dlists[[k]]$location))
    for (i in names(forms[[k]])){
      mml[[k]][[i]] <- model.matrix(forms[[k]][[i]], m[[k]])
      mx[[k]][[i]] <- length(unlist(mx[[k]])) + seq_len(ncol(mml[[k]][[i]][,-1,drop=FALSE]))
    }
    X[[k]] <- compress.model.matrices(mml[[k]])
    npc <- unlist(lapply(mml[[k]], ncol))  - 1
    whichparcov[[k]] <- rep(names(npc), npc)
  }
  if (is.null(pformula)) pformula <- ~1
  Xp <- model.matrix(pformula, m[[1]])[,-1,drop=FALSE]

  ## Convert event internally to numbers based on factor levels
  event <- model.extract(m[[1]], "event")
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
  ntparsl <- nthetal + ncoveffsl # total number of time-to-event related pars for each component
  parindsl <- rep(1:K, ntparsl) # index identifying component for each of those parameters
  ntpars <- sum(ntparsl) # total number of parameters related to time-to-event distributions
  npars <- nppars + ntpars # total number of parameters
  ## identify parameters of particular kinds
  parclass <- rep(c("prob","time"), c(nppars, ntpars))  # mixing prob or time-to-event related
  baseorcov <- c(rep(c("pbase","pcov"), c(length(dists) -1, ncoveffsp)),  # baseline or covariate effect
                 rep(rep(c("tbase","tcov"), K), as.vector(rbind(nthetal, ncoveffsl))))
  parcov <- character(npars) # identify baseline parameter which a covariate effect modifies
  parcov[parclass=="prob"  &  baseorcov=="pcov"] <- rep(paste0("prob",2:K), each=ncovsp)
  parcov[baseorcov=="tcov"] <- unlist(whichparcov)

  ## Build initial values for each parameter type for optimisation
  if (is.null(initp)) initp <- rep(1/K, K)
  alpha <- initp
  theta_inits <- cov_inits <- covparsl <- vector(K, mode="list")
  for (k in 1:K) {
    covparsl[[k]] <- nthetal[k] + seq_len(ncoveffsl[k])
    if (is.null(inits[[k]])) {
      Y <- check.flexsurv.response(model.extract(m[[k]], "response"))
      yy <- ifelse(Y[,"status"]==3 & is.finite(Y[,"time2"]), (Y[,"time1"] + Y[,"time2"])/2, Y[,"time1"])
      yy  <- yy[event==k | is.na(event)]
      # wt <- yy*weights*length(yy)/sum(weights)
      dlists[[k]]$inits <- expand.inits.args(dlists[[k]]$inits)
      ## This extra stuff is needed for the Weibull, to use a  survreg fit to
      ## get "initial values"
      if (dlists[[k]]$name %in% c("exp","weibull.quiet","lnorm","weibullPH","llogis"))
        inits.aux <- list(forms=forms[[k]], data=if(missing(data)) NULL else data, weights=temp$weights,
                          control=sr.control,
                          counting=(attr(model.extract(m[[k]], "response"), "type")=="counting")
        ) else inits.aux <- aux[[k]]
      ## Auto-generate initial values using the heuristic for that distribution
      auto.inits <- dlists[[k]]$inits(t=yy,mf=m[[k]],mml=mml[[k]],aux=inits.aux)
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
  inits_alpha <- if (K==1) numeric() else qmnlogit(inits_probs)
  inits_covp <- rep( rep(0,ncovsp), K-1)   # order by covariate within probability
  names(inits_covp) <- sprintf("prob%s(%s)", rep(2:K, each=ncovsp), rep(colnames(Xp), K-1))
  inits_theta <- numeric()
  for (k in 1:K){
    inits_theta <- c(inits_theta, par.transform(theta_inits[[k]], dlists[[k]]), cov_inits[[k]])
  }

  ## Full loglikelihood function, required for the direct likelihood maximisation
  ## (note this is different from the "complete data" loglikelihood used in EM)
  loglik_flexsurvmix <- function(parsopt, ...){
    pars <- inits_all
    pars[optpars] <- parsopt

    if (K==1)
      probmat <- pmat <- matrix(1, nrow=nobs, ncol=K)
    else {
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
    }
    ## Remaining parameters are time-to-event stuff
    parsl <- split(pars[nppars + seq_len(ntpars)], parindsl)
    theta <- coveffs <- vector(K, mode="list")

    ## Contribution to likelihood for mixing probabilities for those with known events
    llp_event_known <- matrix(0, nrow=nobs, ncol=K)

    # Set membership prob to 0 where there are partially-observed events
    for (i in seq_len(npartials)){
      non_events <- setdiff(1:K, partial_events[[i]])
      probmat[!is.na(partial_event) & partial_event==i, non_events] <- 0
    }

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
                                      list(formula=locform[[k]], data=data, dist=dists[[k]],
                                           anc=anc[[k]], inits=initsk, subset=need_lik,
                                           aux=aux[[k]],
                                           hessian=FALSE,
                                           fixedpars=TRUE))$logliki)
      llp_event_known[,k] <- as.numeric((!is.na(event) & event==k) * log(pmat[,k]))
    }
    logliki <- - (log(rowSums(probmat*liki, na.rm=TRUE)) + rowSums(llp_event_known))
    res <- sum(logliki)
    attr(res, "indiv") <- logliki
    res
  }
  # Parameter dictionary to be completed with the estimates
  distnames <- sapply(dists, function(x){if(is.list(x))x$name else x})
  res <- data.frame(component = c(evnames, rep(evnames[setdiff(1:K, 1)], each=ncovsp),
                                  rep(evnames, ntparsl)),
                    dist = c(rep("",nppars+1), rep(distnames, ntparsl)),
                    terms = c(paste0("prob",1:K),  names(inits_covp),  names(inits_theta)),
                    parclass = c("prob", parclass),
                    baseorcov = c("pbase", baseorcov),
                    parcov = c("", parcov))
  inits_all <- c(inits_alpha, inits_covp, inits_theta)
  if (isTRUE(fixedpars)) fixedpars <- seq_len(npars)
  if (is.logical(fixedpars) && (fixedpars==FALSE)) fixedpars <- NULL
  if (is.character(fixedpars)){
    if (any(!fixedpars %in% res$terms)) {
      bad_fixedpars <- fixedpars[! fixedpars %in% res$terms[-1]]
      stop(paste(bad_fixedpars,collapse=","), " not in model terms")
    }
    fixedpars <- match(fixedpars, as.character(res$terms[-1]))
  }
  if (is.numeric(fixedpars) && any(!fixedpars %in% 1:npars)) 
    stop(sprintf("`fixedpars` should all be in 1,...,%s",npars))
  optpars <- setdiff(seq_len(npars), fixedpars)
  fixed <- (length(optpars) == 0)
  inits_opt <- inits_all[optpars]
  if (is.null(optim.control$fnscale))
    optim.control$fnscale <- 10000

  method <- match.arg(method, c("direct","em"))
  if (!any(is.na(event))) method <- "em" # likelihoods factorise, use EM code with one iteration
  if (method=="direct"){
    if (length(optpars) > 0){
      if (is.null(optim.control$ndeps))
        optim.control$ndeps = rep(1e-06, length(inits_opt))
      opt <- optim(inits_opt, loglik_flexsurvmix, hessian=FALSE, method="BFGS",
                   control=optim.control, ...)
      if (opt$convergence==1)
        warning("Iteration limit in optim() reached without convergence. Reported estimates are not the maximum likelihood. Increase \"maxit\" or simplify the model.")
    } else {
      opt <- list(par=inits_all, value=loglik_flexsurvmix(inits_all))
    }
    logliki <- -attr(loglik_flexsurvmix(opt$par), "indiv")

    ## Transform mixing probs back to natural scale for presentation
    opt_all <- inits_all
    opt_all[optpars] <- opt$par
    res_probs <- if (K>1)  pmnlogit(opt_all[1:(K-1)])  else 1
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
    est.t <- c(NA, opt_all)
    res <- cbind(res, data.frame(
                      est = c(res_probs, res_covp, unlist(res_cpars)),
                      est.t = est.t))

    loglik <- - as.vector(opt$value)
    if (!fixed) {
      # the hessian computation is potentially extremely time consuming!
       cov <- .hess_to_cov(.hessian(loglik_flexsurvmix, opt$par), 
                           hess.control$tol.solve, hess.control$tol.evalues)
    } else {
      cov <- NULL
    }
  }

  else if (method=="em") {
    if (!is.list(em.control)) em.control <- list()
    if (is.null(em.control$reltol)) em.control$reltol <- sqrt(.Machine$double.eps)
    if (is.null(em.control$var.method)) em.control$var.method <- "direct"
    converged <- FALSE
    iter <- 0
    alpha <- inits_alpha
    covp <- inits_covp
    theta <- theta_inits
    covtheta <- cov_inits

    fixedpars_em <- split_fixedpars(fixedpars, nppars, ntparsl, K, nthetal)
    optpars_em <- split_fixedpars(optpars, nppars, ntparsl, K, nthetal)
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
        fs <- flexsurvreg(formula=locform[[k]], data=data, dist=dists[[k]],
                          anc=anc[[k]], inits=theta[[k]], aux=aux[[k]], fixedpars=TRUE,
                          hessian=FALSE)
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
        parsfull <- numeric(nppars)
        parsfull[optpars_em$p] <- pars
        parsfull[fixedpars_em$p] <- c(inits_alpha, inits_covp)[fixedpars_em$p]
        alphamat <- matrix(c(0,parsfull[1:(K-1)]), nrow=nobs, ncol=K, byrow=TRUE)
        if (ncovsp > 0) {
          for (k in 2:K){
            cpinds <- K - 1 + (k-2)*ncovsp + 1:ncovsp
            alphamat[,k] <- alphamat[,k] + Xp %*% parsfull[cpinds]
          }
        }
        pmat <- exp(alphamat)
        pmat <- pmat / rowSums(pmat)
        -sum(w*log(pmat))
      }

      #  if setting a control argument here in  the future, note ndeps has  to be different
      # length from others.
      alpha <- inits_alpha
      covp <- inits_covp
      if (length(fixedpars_em$p) == nppars) { # all prob-related parameters fixed
        hess_full_p <- matrix(nrow=0,ncol=0)
      } else {
        parsp <- optim(c(alpha,covp)[optpars_em$p], loglik_p_em, method="BFGS", hessian=TRUE,
                       control = em.control$optim.p.control)
        if(any(names(optpars_em$p)=="p"))
          alpha[optpars_em$p[names(optpars_em$p)=="p"]] <- parsp$par[names(optpars_em$p)=="p"]  
        if(any(names(optpars_em$p)=="pcov"))
          covp[optpars_em$p[names(optpars_em$p)=="pcov"] - (K-1)] <- parsp$par[names(optpars_em$p)=="pcov"] 
        hess_full_p <- parsp$hessian
      }
      probs <- pmnlogit(alpha)
      
      ## call flexsurvreg for each component on weighted dataset to estimate component-specific pars
      thetanew <- covthetanew <- ttepars <- hess_full_t <- vector(K, mode="list")
      ll <- numeric(K)
    
      ctrl <- optim.control
      for (k in 1:K) {
        if (is.null(optim.control$ndeps))
          ctrl$ndeps = rep(1e-06, length(optpars_em$t[[as.character(k)]]))
        fs <- do.call("flexsurvreg",  # need do.call to avoid environment faff with supplying weights
                      list(formula=locform[[k]], data=data, dist=dists[[k]],
                           anc=anc[[k]], inits=c(theta[[k]],covtheta[[k]]),
                           fixedpars = fixedpars_em$t[[as.character(k)]],
                           weights=w[,k],
                           aux=aux[[k]],
                           subset=w[,k]>0, 
                           control=ctrl))
        ll[k] <- fs$loglik
        thetanew[[k]] <- fs$res[seq_len(nthetal[k]),"est"]
        if (ncoveffsl[k] > 0)
          covthetanew[[k]] <- fs$res[nthetal[k] +  seq_len(ncoveffsl[k]), "est"]
        ttepars[[k]] <- c(thetanew[[k]], covthetanew[[k]])
        hess_full_t[[k]] <- fs$opt$hessian
        if (is.null(hess_full_t[[k]])) hess_full_t[[k]] <- matrix(nrow=0,ncol=0)
      }
      
      theta <- thetanew
      covtheta <- covthetanew
      logliknew <- sum(ll)
      if (!any(is.na(event))) 
        converged <- TRUE
      else if (iter > 0)
        converged <-  (abs(logliknew / loglik - 1) <= em.control$reltol)
      loglik <- logliknew
      est <- c(probs, covp, unlist(ttepars))
      theta.t <- parlist.transform(theta,dlists)
      ttepars.t <- vector(mode="list")
      for (k in 1:K){
        ttepars.t[[k]] <- c(theta.t[[k]], covtheta[[k]])
      }
      est.t <- c(NA, alpha, covp, unlist(ttepars.t))
      if (is.numeric(em.control$trace) && em.control$trace > 0)
        cat(sprintf("loglik=%s\n",loglik))
      if (is.numeric(em.control$trace) && em.control$trace > 1)
        print(est)
      iter <- iter + 1
    }
    res <- cbind(res,
                 data.frame(est=est, est.t=est.t))
    ll <- loglik_flexsurvmix(est.t[-1][optpars])
    loglik <- -as.vector(ll)
    logliki <- -attr(ll, "indiv")
    
    if (em.control$var.method=="direct"){
      hess <- .hessian(loglik_flexsurvmix, est.t[-1][optpars])
      cov <- .hess_to_cov(hess, hess.control$tol.solve, hess.control$tol.evalues)
    }
    else if (em.control$var.method=="louis") 
      cov <- flexsurvmix_louis(K, nthetal, dlists, ncoveffsl, 
                          locform, data, dists, anc,  aux, 
                          fixedpars_em, optpars_em, inits, optim.control, 
                          ttepars.t, nobs, ntparsl, nppars, 
                          w, alpha, covp, ncovsp, Xp,
                          hess_full_p, hess_full_t, hess.control) 
    opt <- NULL
  }

  ## Add standard errors to results data frame, given covariance matrix.
  ## var ( 1 - p1 - p2 - .. ) = var(p1) + var(p2) - cov(p1,p2) - ...
  res$fixed <- c(NA, rep(FALSE, npars))
  res$fixed[1 + fixedpars] <- TRUE
  optp <- optpars[optpars < K]
  res$se <- rep(NA, npars + 1)
  if (!fixed){
      if (length(optp) > 0){
        covp <- cov[optp, optp, drop=FALSE]
        res$se[1] <- if (K==1) NA else sqrt(sum(diag(covp)) - sum(covp[lower.tri(covp)]))
      }
      res$se[1 + optpars] <- sqrt(diag(cov))
  } 
  names.first <- c("component","dist","terms","est","est.t","se")
  res <- res[,c(names.first, setdiff(names(res), names.first)),drop=FALSE]
  rownames(res) <- NULL

  mcomb <- do.call("cbind", m)
  mcomb <- mcomb[,!duplicated(names(mcomb))]
  attr(mcomb, "covnames") <- unique(unlist(lapply(m, function(x)attr(terms(x),"term.labels"))))
  attr(mcomb, "covnames.main") <- unique(unlist(lapply(m, function(x)rownames(attr(terms(x),"factors"))[-1])))
  
  res <- list(call=match.call(),
              res=res, loglik=loglik,  cov=cov,
              npars=npars, AIC=-2*loglik + 2*npars,
              K=K, dists=distnames,  dlists=dlists, dfns=dfns,
              opt=opt,
              fixedpars=fixedpars, optpars=optpars,
              logliki=logliki,
              evnames=evnames,
              nthetal=nthetal, parindsl=parindsl,
              ncoveffsl=ncoveffsl,ncoveffsp=ncoveffsp, covparsl=covparsl,
              all.formulae=forms, pformula=pformula,
              data=list(mf=m, mfcomb=mcomb), mx=mx)
  class(res) <- "flexsurvmix"
  res
}

check.formula.flexsurvmix <- function(formula, dlist){
  if (!is.list(formula)){
    if (!inherits(formula,"formula")) stop("\"formula\" must be a formula object")
    formula <- list(formula)
  }
  for (i in seq_along(formula)){
    if (!inherits(formula[[i]],"formula"))
      stop(sprintf("\"formula[[%s]]\" must be a formula object", i))
    labs <- attr(terms(formula[[i]]), "term.labels")
    if (!("strata" %in% dlist$pars)){
      strat <- grep("strata\\((.+)\\)",labs)
      if (length(strat) > 0){
        cov <- gsub("strata\\((.+)\\)","\\1",labs[strat[1]])
        warning("Ignoring \"strata\" function: interpreting \"",cov, "\" as a covariate on \"", dlist$location, "\"")
      }
    }
    if (!("frailty" %in% dlist$pars)){
      fra <- grep("frailty\\((.+)\\)",labs)
      if (length(fra) > 0){
        warning("frailty models are not supported and behaviour of frailty() is undefined")
      }
    }
  }
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
        stop(sprintf("names(%s) = %s, but this should match levels(factor(event)) = %s",
                argname, paste(narg, collapse=","), paste(evnames,collapse=",")))
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
    keep.cols <- c("component","dist","terms","est","est.t","se")
    res <- x$res[,keep.cols,drop=FALSE]
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
  probs <- if (K==1) 1 else pmnlogit(x$res$est.t[2:K])
  pcov <- x$res$est.t[grep("prob[[:digit:]]+\\(.+\\)", x$res$terms)]
  c(probs, pcov, unlist(est))
}

##' @export
model.frame.flexsurvmix <- function(formula, ...)
{
  x <- formula
  x$data$m
}

# Multinomial logistic transform and inverse transform

# returns transformed probs relative to first component, excluding first component
qmnlogit <- function(probs){
  if  (length(probs) == 1)
    qlogis(probs)
  else
    log(probs[-1] / probs[1])
}

# returns probs including first component, given transformed probs with first component excluded
pmnlogit <- function(alpha){
  c(1, exp(alpha)) / (1 + sum(exp(alpha)))
}

## Convert "fixedpars" (vector of parameter indices to fix) from 
## joint full MLE form (indices into full parameter vector) to EM form 
## (indices of parameters in each submodel to be maximised in the M step)
## e.g. 
#nppars <- 4    # number of prob pars (including cov effects)
#ntparsl <- c(3, 3, 4) # number of time to event pars for each event (including cov effects)
#indices  1,2,3,4  are prob pars, 5,6,7 time to event 1 pars
# 8,9,10 are time to event 2 pars, and 11,12,13,14 are time to event 3 pars
#fixedpars <- c(2,3,   6,7,    9,   12)
#correspond to indices 2,3    2,3     2   2)  in the four submodels 
# so function returns list(p=c(2,3), t=list(c(2,3), 2, 2))
# List component will be numeric(0) for prob model, or NULL for time model, if no pars fixed for that model
# names of time models will be "1","2",.

split_fixedpars <- function(fixedpars, nppars, ntparsl, K, nthetal){
  fixedpars_p <- fixedpars[fixedpars <= nppars]
  if (!is.null(fixedpars_p))
    names(fixedpars_p) <- ifelse(fixedpars_p <= K-1, "p", "pcov")
  ft <- fixedpars[fixedpars >= nppars] - nppars
  fixedpars_t <- 
    split(sequence(ntparsl)[ft],  
          rep(seq_along(ntparsl), ntparsl)[ft])
  for (i in 1:K){
    ic <- as.character(i)
    if (!is.null(fixedpars_t[[ic]]))
      names(fixedpars_t[[ic]]) <- ifelse(fixedpars_t[[ic]] <= nthetal[i], "theta",  "covtheta")
    else fixedpars_t[[ic]] <- numeric()
  }
  fixedpars_t <- fixedpars_t[as.character(1:K)]
  list(p=fixedpars_p, t=fixedpars_t)
}
