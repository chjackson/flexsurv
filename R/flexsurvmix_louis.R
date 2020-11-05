## Calculate the covariance matrix for a flexsurvmix model fitted by the EM algorithm
## using the method of Louis (JRSSB, 1982)

flexsurvmix_louis <- function(K, nthetal, dlists, ncoveffsl, 
                             locform, data, dists, anc,  aux, 
                             fixedpars_em,  optpars_em, inits, optim.control, 
                             ttepars.t, nobs, ntparsl, nppars, 
                             w, alpha, covp, ncovsp, Xp, 
                             hess_full_p, hess_full_t, hess.control) 
{ 
  ESST_t <- vector(K, mode="list")
  ctrl <- optim.control
  for (k in 1:K){ 
    ## unweighted indiv loglik (not minus loglik)
    if (is.null(optim.control$ndeps)) ctrl$ndeps = rep(1e-06, length(optpars_em$t[[as.character(k)]]))
    logliki_flexsurvreg <- function(parsopt){
      pars <- inits[[k]]
      pars[optpars_em$t[[k]]] <- parsopt
      theta <- inv.transform(pars[seq_len(nthetal[k])], dlists[[k]])
      coveffs <- pars[nthetal[k] + seq_len(ncoveffsl[k])]
      initsk <- c(theta, coveffs)
      fs <- do.call("flexsurvreg", list(formula=locform[[k]], 
                                        data=data, 
                                        dist=dists[[k]],
                                        anc=anc[[k]], 
                                        inits=initsk, 
                                        aux=aux[[k]], 
                                        fixedpars = fixedpars_em$t[[as.character(k)]],
                                        hessian=FALSE, 
                                        control=ctrl)
                    )
      fs$logliki
    }
    parsopt <- ttepars.t[[k]][optpars_em$t[[k]]]
    nopt <- length(parsopt)
    J <- numDeriv::jacobian(logliki_flexsurvreg, parsopt) # nindiv x npars. 
    ## Quick with nindiv=1000
    SST <- array(dim=c(nobs, nopt, nopt))
    for (i in 1:nopt) { 
      for (j in 1:nopt) {
        SST[,i,j]  <- J[,i] * J[,j]
      }
    }
    # w[,k] # nindiv vector 
    ESST_t[[k]] <- apply(SST * w[,k], c(2,3), sum)
  }

  # Analytic derivatives of imputed-data log-likelihood for the component 
  # membership probs (and covariates on these)
  alphamat <- matrix(c(0,alpha), nrow=nobs, ncol=K, byrow=TRUE)  # by individual
  if (ncovsp > 0) {
    for (k in 2:K){
      cpinds <-  (k-2)*ncovsp + 1:ncovsp
      alphamat[,k] <- alpha[k-1] + Xp %*% covp[cpinds]
    }
  }
  dmat <- array(dim=c(nobs, nppars, K))
  sumexp <- rowSums(exp(alphamat))  
  for (j in 2:K) {
    esum <- exp(alphamat[,j]) / sumexp
    for (k in 1:K) { 
      dmat[,j-1,k] <- if (k==1) - esum else if (k==j) 1 - esum else (-exp(alphamat[,j]) / sumexp)
      if (ncovsp > 0)  {
        for (i in 1:ncovsp){
          ind <- K-1 + (j-2)*ncovsp + i
          dmat[,ind,k] <- if (k==1) - Xp[,i] * esum else if (k==j) Xp[,i] * (1 - esum) else Xp[,i] * (-exp(alphamat[,j]) / sumexp)
        }
      }
    }
  }
  dmat <- -dmat
  Dw <- array(dim=c(nobs, nppars, K))
  DDw <- array(dim=c(nobs, nppars, nppars, K))
  for (k in 1:K){
    for (i in 1:nppars){
      Dw[,i,k] <- dmat[,i,k] * w[,k]
      for (j in 1:nppars){
        DDw[,i,j,k] <- dmat[,i,k] * dmat[,j,k] * w[,k]
      }
    }
  }
  SDw <- array(dim=c(nobs, nppars))
  for (i in 1:nppars){
    SDw[,i] <- rowSums(Dw[,i,])
  }
  SDw <- SDw[,optpars_em$p,drop=FALSE]
  DDw <- DDw[,optpars_em$p,optpars_em$p,,drop=FALSE]
  noptp <- length(optpars_em$p)
  ESST_p <- array(dim=c(noptp, noptp))
  for (i in 1:noptp){
    for (j in 1:noptp){
      ESST_p[i,j] <- sum(DDw[,i,j,]) + sum(SDw[,i])*sum(SDw[,j]) - sum(SDw[,i]*SDw[,j])
    }
  }
  
  ESST <- Matrix::bdiag(ESST_p, Matrix::bdiag(ESST_t))
  hess_full <- - Matrix::bdiag(hess_full_p, Matrix::bdiag(hess_full_t))
  hess <- hess_full - ESST
  cov <- .hess_to_cov(-hess, hess.control$tol.solve, hess.control$tol.evalues)
  cov # of loglik
}

