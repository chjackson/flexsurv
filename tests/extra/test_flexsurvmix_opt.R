
if (em_supported){
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"),  method="em")
  xd <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                    dists=c("gamma","gamma"),  method="direct")
  
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"), pformula = ~x+y,  method="em")
  xd <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                    dists=c("gamma","gamma"), pformula = ~x+y,  method="direct")
  
  x <- flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                   dists=c("gamma","gamma"),  method="em", fixedpars=1)
  
}