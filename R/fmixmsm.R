##' Constructor for a mixture multi-state model based on flexsurvmix
##'
##' @param ... Named arguments.  Each argument should be a fitted model as
##'   returned by \code{\link{flexsurvmix}}.  The name of each argument names
##'   the starting state for that model.
##'
##' @return A list of \code{\link{flexsurvmix}} objects, with the following
##'   attribute(s):
##'
##'   \code{pathways} A list of all potential pathways until absorption, for
##'   models without cycles.   For models with cycles this will have an element
##'   \code{has_cycle=TRUE}, plus the pathways discovered before the function
##'   found the cycle and gave up.
##'
##' @export
fmixmsm <- function(...){
  args <- list(...)
  starts <- names(args) 
  evlist <- lapply(args, function(x)x$evnames)
  names(evlist) <- starts
  plist <- list(transient_unvisited =  names(evlist),
              transient_visited = character(),
              pathways = list(),
              pathway_current = character())
  pathways <- get_pathways(evlist[1], evlist, plist)
  ret <- args
  attr(ret, "pathways") <- pathways$pathways
  attr(ret, "pathway_str") <- sapply(attr(ret, "pathways"), function(x)paste(x,collapse="-"))
  ret
}

get_pathways <- function(mod_current, mods, ret){
  ## TODO error handling 
  if (isTRUE(ret$has_cycle)) return(ret)
  fromstate <- names(mod_current)
  tostates <- mod_current[[1]]
  absorbing <- setdiff(unlist(mods),names(mods))
  ret$pathway_current <- c(ret$pathway_current, fromstate)
  ret$transient_visited <- c(ret$transient_visited, fromstate)
  ret$transient_unvisited <- setdiff(ret$transient_unvisited, fromstate)
  nd <- length(tostates)
  pcurr <- ret$pathway_current 
  for (j in 1:nd){
    if (tostates[j] %in% ret$transient_visited){
      return(list(has_cycle=TRUE))  # could we carry on, ignore cycles, return all paths to absorption?  don't do unless we need. 
    } 
    else if (tostates[j] %in% ret$transient_unvisited){
      ret <- get_pathways(mods[tostates[j]], mods, ret)
    }
    else if (tostates[j] %in% absorbing){
      ret$pathways <- c(ret$pathways, list(c(pcurr, tostates[j])))
    } else stop("Shouldn't reach here, please report a bug")
  }
  ret
}

##' Probability of each pathway taken through a mixture multi-state model
##'
##'
##' @param x Object returned by \code{\link{fmixmsm}}, representing a multi-state
##'   model built from piecing together mixture models fitted by
##'   \code{\link{flexsurvmix}}.
##'
##' @param final If \code{TRUE} then the probabilities of pathways with the same
##'   final state are added together, to produce the probability  of each
##'   ultimate outcome or absorbing state from the multi-state model.
##'
##' @inheritParams mean_flexsurvmix
##' 
##' @return Data frame of pathway probabilities by covariate value and pathway.
##' 
##'
##' @export
prob_pathway <- function(x, newdata=NULL, final=FALSE, B=NULL){
    pathways <- attr(x, "pathways")
    if (isTRUE(pathways$has_cycle))
        stop("models with cycles not supported in this function")
    nmods <- length(x)
    probs <- vector(nmods, mode="list")
    names(probs) <- names(x)
    for (i in seq_along(x)){
        probs[[i]] <- probs_flexsurvmix(x[[i]], newdata=newdata)
    }
    npaths <- length(pathways)  
    ppath <- vector(npaths, mode="list")
    for (p in seq_along(pathways)){
        plen <- length(pathways[[p]])
        ppath[[p]] <- 1
        for (i in 1:(plen-1)){
          cur_state <- pathways[[p]][i]
          next_state <- pathways[[p]][i+1]
          cur_prob <- probs[[cur_state]] 
          ppath[[p]] <- ppath[[p]] * cur_prob$val[cur_prob$event==next_state] 
        }
    }
    finalstate <- sapply(attr(x, "pathways"), function(x)x[length(x)])
    ncovs <- if(is.null(newdata)) 1 else nrow(newdata)
    finalstate <- rep(finalstate, each=ncovs)
    pathway <- rep(attr(x, "pathway_str"), each=ncovs)
    ppath <- data.frame(final=finalstate,
                        pathway=pathway, 
                        val = unlist(ppath))
    if (!is.null(newdata)) { 
      nd <- newdata[rep(seq_len(ncovs), npaths),]
      ppath <- cbind(nd, ppath)
    }
    rownames(ppath) <- NULL
    if (final) { 
      ppath <- ppath %>% 
        dplyr::group_by_at(c(names(newdata), "final")) %>% 
        dplyr::summarise(val=sum(.data$val))
    } 
    
    if (is.numeric(B) && B > 1){
      res <- matrix(nrow=B, ncol=length(ppath$val))
      res[1,] <- ppath$val
      for (i in 2:B){
        xrep <- resample_pars_fmixmsm(x)
        res[i,] <- prob_pathway(xrep, newdata=newdata, final=final, B=NULL)$val
      }
      resci <- apply(res, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
      ppath$lower <- resci[1,]
      ppath$upper <- resci[2,]
    }
    as.data.frame(ppath)
}

resample_pars_fmixmsm <- function(x){
  xrep <- lapply(x, resample_pars)
  attributes(xrep) <- attributes(x)
  xrep
}

##' Mean time to final state in a mixture multi-state model
##'
##' Calculate the mean time from the start of the process to a final (or
##' "absorbing") state in a mixture multi-state model.  Models with cycles are
##' not supported.
##'
##' @inheritParams prob_pathway
##'
##' @param final If \code{TRUE} then the mean time to the final state is
##'   calculated for each final state, by taking a weighted average of the mean
##'   time to travel each pathway ending in that final state, weighted by the
##'   probability of the pathway.   If \code{FALSE}  (the default) then a
##'   separate mean is calculated for each pathway.
##'
##' @return A data frame of mean times to absorption, by covariate values and
##'   pathway (or by final state)
##' 
##' @export
mean_tofinal <- function(x, newdata=NULL, final=FALSE, B=NULL){
  pathways <- attr(x, "pathways")
  nmods <- length(x)
  means <- vector(nmods, mode="list")
  names(means) <- names(x)
  for (i in seq_along(x)){
    means[[i]] <- mean_flexsurvmix(x[[i]], newdata=newdata)
  }
  npaths <- length(pathways)  
  meanp <- vector(npaths, mode="list")
  for (p in seq_along(pathways)){
    plen <- length(pathways[[p]])
    meanp[[p]] <- 0
    for (i in 1:(plen-1)){
      cur_state <- pathways[[p]][i]
      next_state <- pathways[[p]][i+1]
      cur_mean <- means[[cur_state]]
      meanp[[p]] <- meanp[[p]] + cur_mean$val[cur_mean$event==next_state] 
    }
  }    
  finalstate <- sapply(attr(x, "pathways"), function(x)x[length(x)])
  ncovs <- if(is.null(newdata)) 1 else nrow(newdata)
  finalstate <- rep(finalstate, each=ncovs)
  pathway <- rep(attr(x, "pathway_str"), each=ncovs)
  meanp <- data.frame(final=finalstate,
                      pathway=pathway, 
                      val = unlist(meanp))
  if (!is.null(newdata)) { 
    nd <- newdata[rep(seq_len(ncovs), npaths),]
    meanp <- cbind(nd, meanp)
  }
  rownames(meanp) <- NULL
  if (final) {  
    probs <- prob_pathway(x=x, newdata=newdata, final=FALSE, B=NULL) %>% dplyr::rename(prob="val")
    meanp <- meanp %>% 
      dplyr::left_join(probs, by=c(names(newdata), "pathway", "final")) %>%
      dplyr::group_by_at(c(names(newdata), "final")) %>% 
      dplyr::summarise(val=sum(.data$val*.data$prob)/sum(.data$prob))
  } 
  if (is.numeric(B) && B > 1){
    res <- matrix(nrow=B, ncol=length(meanp$val))
    res[1,] <- meanp$val
    for (i in 2:B){
      xrep <- resample_pars_fmixmsm(x)
      res[i,] <- mean_tofinal(xrep, newdata=newdata, final=final, B=NULL)$val
    }
    resci <- apply(res, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
    meanp$lower <- resci[1,]
    meanp$upper <- resci[2,]
  }
  as.data.frame(meanp)
}


##' Quantiles of the distribution of the time until reaching a final state in a mixture multi-state model
##' 
##' 
##' @inheritParams prob_pathway
##' 
##' @param final If \code{TRUE} then the mean time to the final state is
##'   calculated for each final state, by taking a weighted average of the mean
##'   time to travel each pathway ending in that final state, weighted by the
##'   probability of the pathway.   If \code{FALSE}  (the default) then a
##'   separate mean is calculated for each pathway.
##'
##' @param n Number of individual-level simulations to use to characterise the
##'   time-to-event distributions
##'
##' @param probs Quantiles to calculate, by default, \code{c(0.025, 0.5, 0.975)}
##'
##' @return Data frame of quantiles of the time to final state by pathway and
##'   covariate value, or by final state and covariate value.
##'   
quantile_tofinal <- function(x, newdata=NULL, final=FALSE, B=NULL, n=10000, probs=c(0.025, 0.5, 0.975)){
  pathways <- attr(x, "pathways")
  nmods <- length(x)
  sims <- vector(nmods, mode="list")
  names(sims) <- names(x)
  if (final) {  
    ppath <- prob_pathway(x=x, newdata=newdata, final=FALSE, B=NULL) %>% 
      dplyr::rename(prob="val") %>%
      dplyr::mutate(n = round(n*prob))
  }
  for (i in seq_along(x)){
    sims[[i]] <- simt_flexsurvmix(x[[i]], newdata=newdata, n=n)
  }
  npaths <- length(pathways)  
  simsum <- vector(npaths, mode="list")
  for (p in seq_along(pathways)){
    plen <- length(pathways[[p]])
    sm <- 0
    for (i in 1:(plen-1)){
      cur_state <- pathways[[p]][i]
      next_state <- pathways[[p]][i+1]
      sm <- sm + sims[[cur_state]][,next_state]
    }
    simsumdf <- sims[[1]][,colnames(newdata)]
    simsumdf$pathway <- attr(x,"pathway_str")[[p]]
    simsumdf$sm <- sm
    
    ## TESTME
    if (final){
      simsum[[p]] <- simsumdf %>%
        dplyr::left_join(ppath, by=c(colnames(newdata))) %>%
        dplyr::group_by_at(c(colnames(newdata))) %>%
        ## Keep only the first n of the sampled rows
        ## where n is weighted by the prob of the pathway
        dplyr::group_modify(~{.x[1:.x$n[1],]}) 
    } else {
      simsum[[p]] <- simsumdf %>%
      dplyr::group_by_at(c("pathway", colnames(newdata)))  %>%
      dplyr::summarise(probs=probs,
                       val = quantile(.data$sm, p=probs,na.rm=TRUE))
    }
  }
  resq <- do.call("rbind", simsum)
  
  if (final){
    resq <- resq %>%
      dplyr::group_by_at(c("final", colnames(newdata)))  %>%
      dplyr::summarise(probs=probs,
                       val = quantile(.data$sm, p=probs,na.rm=TRUE))
  }

  rownames(resq) <- NULL
  
  if (is.numeric(B) && B > 1){
    res <- matrix(nrow=B, ncol=length(resq$val))
    res[1,] <- resq$val
    for (i in 2:B){
      xrep <- resample_pars_fmixmsm(x)
      res[i,] <- quantile_tofinal(xrep, newdata=newdata, final=final, B=NULL, n=n, probs=probs)$val
    }
    resci <- apply(res, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
    resq$lower <- resci[1,]
    resq$upper <- resci[2,]
  }
  as.data.frame(resq)
}

##' Simulate times to competing events from a mixture multi-state model
##'
##' @inheritParams probs_flexsurvmix
##'
##' @param n Number of simulations
##'
##' @return Data frame with \code{n*m} rows and a column for each competing
##'   event, where \code{m} is the number of alternative covariate values, that
##'   is the number of rows of \code{newdata}.   The simulated time represents
##'   the time to that event conditionally on that event being the one that
##'   occurs.  This function doesn't simulate which event occurs.
##'   
simt_flexsurvmix <- function(x, newdata=NULL, n){
  if (is.null(newdata)) newdata <- default_newdata(x)
  if (is.null(newdata)) {
    ncovvals <- 1
  }
  else {
    newdata <- as.data.frame(newdata)
    ncovvals <- nrow(newdata)
  }
  simdf <- as.data.frame(matrix(nrow=n*ncovvals, ncol=x$K))
  names(simdf) <- x$evnames
  for (k in 1:x$K){
    pars <- get_basepars(x, newdata, k)
    pars$n <- n*ncovvals
    simdf[[k]] <- do.call(x$dfns[[k]]$r, pars)
  }
  if (!is.null(newdata)){
    newdatarep <- newdata[rep(seq_len(ncovvals),n),,drop=FALSE]
  }
  simdf <- cbind(newdatarep, simdf)
  simdf <- simdf[do.call("order", newdatarep),,drop=FALSE]
  simdf
}
