
if (interactive()){
  fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
  cs <- coxsnell_flexsurvreg(fitg)
  
  ## Model doesn't appear to fit well since the cumulative hazards are underestimated. 
  ## In this example, this is because the dataset is small, hence the point estimate is noisy. 
  plot(cs$"(qexp)", cs$est, pch=19, xlab="Theoretical quantiles", ylab="Cumulative hazard")
  abline(a=0,b=1,col="red",lwd=2)
  
  ## Alternative way to produce the same plot using "qqplot" 
  qy <- qexp(ppoints(nrow(cs),0))
  qqplot(qy, cs$est)
  abline(a=0,b=1, col="red", lwd=2)

  ## A log transform may or may not bring out the pattern more clearly
  plot(log(cs$"(qexp)"), log(cs$est), pch=19)
  abline(a=0,b=1, col="red", lwd=2)
  
  ## In the model `fitg`, the fitted cumulative hazard is lower than the true cumulative hazard 
  ## Another way to show this
  plot(fitg, type="cumhaz", ci=FALSE)

  ## Alternative situation where the true model is fitted to simulated data 
  y <- rweibull(10000, 2, 2)
  fite <- flexsurvreg(Surv(y) ~ 1, dist="weibull")
  cs <- coxsnell_flexsurvreg(fite)
  ## The model fits well
  plot(cs$"(qexp)", cs$est, pch=19, xlab="Theoretical quantiles", ylab="Cumulative hazard")
  abline(a=0,b=1,col="red",lwd=2)
}