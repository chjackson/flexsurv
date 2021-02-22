
test_that("RMST/Mean/Median calculations are working",{
  
  fs1 = flexsurvreg(Surv(rectime, censrec)~group ,dist="weibull",data=bc)
  fs2 = flexsurvreg(Surv(rectime, censrec)~group ,dist="exp",data=bc)
  fs3 = flexsurvreg(Surv(rectime, censrec)~group ,dist="llogis",data=bc)
  fs4 = flexsurvreg(Surv(rectime, censrec)~group ,dist="lnorm",data=bc)
  suppressWarnings({  # warnings from temporary overflow during optimisation 
    fs5 = flexsurvreg(Surv(rectime, censrec)~group ,dist="gamma",data=bc)
  })
  res1 = summary(fs1,t=c(Inf),start=0,type="rmst")
  res2 = summary(fs1,type="mean")
  
  res1_len = length(res1)
  for(i in seq_len(res1_len)){
    expect_equal(
      res1[[i]]$est,
      res2[[i]]$est,
      tolerance=1e-3
    )
  }
  
  # Exponential analytical RMST should be consistent w/ analytical
  # mean.
  expect_equal(summary(fs2,type="mean",tidy=T)$est,summary(fs2,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  
  # Analytical mean should closely match result from integration
  expect_equal(summary(fs1,type="mean",tidy=T)$est,summary(fs1,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  expect_equal(summary(fs3,type="mean",tidy=T)$est,summary(fs3,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  expect_equal(summary(fs4,type="mean",tidy=T)$est,summary(fs4,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  expect_equal(summary(fs5,type="mean",tidy=T)$est,summary(fs5,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  
  
  # RMST of exponential to 100 starting at 0 should be the same
  # as RMST to 200 starting at 100.
  res3 = summary(fs2,t=c(100,200),start=c(0,100),type="rmst")
  res3_len = length(res3)
  for(i in seq_len(res1_len)){
    expect_equal(
      res3[[i]]$est[1],
      res3[[i]]$est[2],
      tolerance=1e-3
    )
  }
  
  expect_warning(
    summary(fs1,t=10,type="mean"),
    "Mean selected, but time specified.  For restricted mean, set type to 'rmst'."
  )
  
  expect_warning(
    summary(fs1,t=10,type="median"),
    "Median selected, but time specified."
  )
  
})
