standsurv.flexsurvreg <- function(object, newdata = NULL, at = list(), atreference = 1, type = "survival", t = NULL,
                                  ci = FALSE, B = 1000, cl =0.95, contrast = NULL, seed = NULL) {
  x <- object
  
  if(!is.null(seed)) set.seed(seed)
  
  ## Add checks
  ## Currently restricted to survival or rmst 
  type <- match.arg(type, c("survival", "rmst"))
  contrast <- match.arg(contrast, c("difference", "ratio"))
  ## Currently does not calculate CIs 
  
  ## Check that at is a list and that all elements of at are lists
  if(!is.list(at)){
    stop("'at' must be at list")
  }
  if(any(!sapply(at, is.list))){
    stop("All elements of 'at' must be lists")
  }
  
  ## Contrast numbers
  cnums <- (1:length(at))[-atreference]
  
  ## Standardise over fitted dataset by default
  if(is.null(newdata)){
    data <- model.frame(x)
  } else{
    data <- newdata
  }

  ## If at is not specified then no further manipulation of data is required, we standardise over original or passes dataset
  stand.pred.list <- list()
  for(i in 1:length(at)){
    dat <- data
    covs <- at[[i]]
    covnames <- names(covs)
    for (j in 1:length(covnames)) dat[, covnames[j]] <- covs[j]
 
    pred <- summary(object, type = type, tidy = T, newdata=dat, t=t, ci=F) ## this gives predictions based on MLEs (no bootstrapping for point estimates)
    predsum <- pred %>% group_by(time) %>% summarise("at{i}" := mean(est))
    
    if(ci == TRUE){
  
      if(i==1)      rawsim <- normboot.flexsurvreg(object, B=B, raw=T) ## only run this once, not for every specified _at
  
      X <- form.model.matrix(object, as.data.frame(dat), na.action=na.pass)
      sim.pred <- normbootfn.flexsurvreg(object, t=t, start=0, X=X, fn=summary.fns(object, type), B=B, rawsim=rawsim) # pts, sims, times
      
      stand.pred <- apply(sim.pred, c(2,3), mean)
      stand.pred.quant <- apply(stand.pred, 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE) )
      stand.pred.quant <- as_tibble(t(stand.pred.quant)) %>% rename("at{i}_lci" := "2.5%", "at{i}_uci" := "97.5%")
      
      predsum <- predsum %>% bind_cols(stand.pred.quant)
      stand.pred.list[[i]] <- stand.pred
    }
    
    if(i == 1) {
      standpred <- predsum
    } else {
      standpred <- standpred %>% inner_join(predsum, by ="time")
    }
    
  }
  
  if(contrast == "difference"){
    for(i in cnums){
      standpred <- standpred %>% mutate("contrast{i}_{atreference}" := .data[[paste0("at", i)]] - .data[[paste0("at", atreference)]])
      if(ci == TRUE){
        stand.pred.quant <- apply(stand.pred.list[[i]] - stand.pred.list[[atreference]], 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        stand.pred.quant <- as_tibble(t(stand.pred.quant)) %>% rename("contrast{i}_{atreference}_lci" := "2.5%", "contrast{i}_{atreference}_uci" := "97.5%")
        standpred <- standpred %>% bind_cols(stand.pred.quant)
      }
    }
  }
  if(contrast == "ratio"){
    for(i in cnums){
      standpred <- standpred %>% mutate("contrast{i}_{atreference}" := .data[[paste0("at", i)]] / .data[[paste0("at", atreference)]])
      if(ci == TRUE){
        stand.pred.quant <- apply(stand.pred.list[[i]] / stand.pred.list[[atreference]], 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        stand.pred.quant <- as_tibble(t(stand.pred.quant)) %>% rename("contrast{i}_{atreference}_lci" := "2.5%", "contrast{i}_{atreference}_uci" := "97.5%")
        standpred <- standpred %>% bind_cols(stand.pred.quant)
      }
    }
  }
  
  standpred
}


