# Code from survival package.  Credit to Terry Therneau.
untangle.specials <- function (tt, special, order = 1){
  spc <- attr(tt, "specials")[[special]]
  if (length(spc) == 0) 
    return(list(vars = character(0), terms = numeric(0)))
  facs <- attr(tt, "factors")
  fname <- dimnames(facs)
  ff <- apply(facs[spc, , drop = FALSE], 2, sum)
  list(
    vars = (fname[[1]])[spc],
    terms = seq(ff)[ff & match(attr(tt,"order"), order, nomatch = 0)]
  )
}

# Extract string containing terms for a given ancillary parameter
parseFormulas <- function(formula, locname, ancnames, anc = NULL, data = NULL){
  
  # Capture formula environment
  env <- environment(formula)
  
  # Extract response
  response <- as.character(formula)[2]
  
  # Name vectors of special/parameter names
  names(ancnames) <- ancnames
  specials <- c(locname, ancnames, "strata")
  names(specials) <- c(locname, ancnames, "strata")
  
  # Extract terms
  terms <- terms(formula, special = specials)
  specialTerms <- lapply(specials, untangle.specials, tt = terms)
  
  # Extract location terms
  specialIndices <- do.call(c, lapply(specialTerms, function(x) x$terms))
  if(length(specialIndices) == 0) locTermsPrimary <- attributes(terms)$term.labels
  else locTermsPrimary <- attributes(terms[-specialIndices])$term.labels
  
  # Extract strata terms
  strataTerms <- specialTerms$strata$vars
  pattern <- paste0("strata\\((.+)\\)")
  
  strataVars <- unlist(strsplit(
    gsub(pattern, "\\1", grep(pattern,strataTerms,value=TRUE)),
    ",",
    fixed = TRUE
  ))
  
  
  selector <- as.logical(lapply(
    strataVars,
    function(x){
      val <- eval(parse(text=x), envir = data, enclos = env)
      is.numeric(val)
    }
  ))
  if(any(selector)) warning("Ignoring numeric variables in strata")
  strataVars <- strataVars[!selector]
  if(length(strataVars) == 0) strataTerms <- NULL
  else strataTerms <- paste(strataVars[!selector], collapse = "*")
  
  
  # Extract covariate terms
  covTerms <- lapply(
    specials[-length(specials)],
    function(x){
      labs <- specialTerms[[x]]$vars
      pattern <- paste0(x,"\\((.+)\\)")
      labs <- grep(pattern,labs,value=TRUE)
      if (length(labs)==0) return(NULL)
      labs <- gsub(pattern, "\\1", labs)
      labs
    }
  )
  
  # Construct ancillary formulas
  if(is.null(anc)){
    # If anc is not supplied, construct from formula
    ancTerms <- lapply(
      ancnames,
      function(x){
        ancTerms <- append(covTerms[[x]], strataTerms)
        if(length(ancTerms) == 0) NULL
        else reformulate(ancTerms)
      }
    )
  }else{
    # Otherwise, use anc and warn user if he has also placed ancillary in formula
    lapply(
      ancnames,
      function(x){
        ancTerms <- covTerms[[x]]
        # Throw error if user provides anc and also has ancillary covariates in main formula
        if(length(ancTerms) != 0)
          stop("Covariates for ancillary formulas may not be provided in main formula when 'anc' arugment provided.")
      }
    )
    # Throw error if user tries to put location formula in anc.
    if(!is.null(anc[[locname]]))
      stop("Formula for location parameter may not be provided in 'anc' argument.")
    ancTerms <- lapply(
      ancnames,
      function(x){
        if(!is.null(anc[[x]])){
          ancTerms <- append(attributes(terms(anc[[x]]))$term.labels, strataTerms)
          if(length(ancTerms) == 0) return(NULL)
          else return(reformulate(ancTerms))
        }
        else return(NULL)
      }
    )
  }
  
  # Construct location formula
  locTerms <- append(append(locTermsPrimary, covTerms[[locname]]), strataTerms)
  if(length(locTerms) == 0) locTerms <- "1"
  locForm <- list(as.formula(paste(response, "~", paste(locTerms, collapse="+"))))
  names(locForm) <- locname
  
  # Combine into single list and restore environment
  formList <- lapply(append(locForm, ancTerms), function(x) {
    if(!is.null(x)) environment(x) <- env
    x
  })
  # Remove NULL elements
  formList <- formList[!sapply(formList, is.null)]
  
  return(formList)
}