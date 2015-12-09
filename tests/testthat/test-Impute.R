context("Impute")


test_that("Impute_creation_and_extract_sim",{
  #also look at test-SimFit.R
  set.seed(202)
  
  sim <- SimulateComplete(study.time=365,number.subjects=100,
                          event.rates=c(0.01,0.005),dispersions=0.25)
  
  sim.dropout <- SimulateDropout(sim,drop.mechanism=ConstantRateDrop(rate=0.0025)) 
  fit <- Simfit(sim.dropout,equal.dispersion=TRUE)
  
  impute <- Impute(fit,impute.mechanism = weighted_j2r(trt.weight=0),N=10)
  
  expect_equal("ImputeSim",class(impute))
  expect_equal(summary(sim.dropout)$number.dropouts,impute$dropout)
  expect_equal(2,nrow(impute$imputed.values))
  expect_equal(10,ncol(impute$imputed.values))
  expect_equal(impute$impute.mechanism,weighted_j2r(0))
  
  expect_error(GetImputedDataSet(imputeSim=sim.dropout,index=1))
  expect_error(GetImputedDataSet(imputeSim=impute,index=0))
  expect_error(GetImputedDataSet(imputeSim=impute,index=11))
  expect_error(GetImputedDataSet(imputeSim=impute,index=4.5))
  
  new.sim <- GetImputedDataSet(impute,index=4)
  expect_equal("SingleSim",class(new.sim))
  expect_equal("imputed",new.sim$status)
  
  expect_equal(weighted_j2r(0),new.sim$impute.mechanism)
  expect_equal(new.sim$data$actual.censored.time,sim.dropout$data$censored.time)
  expect_equal(rep(365,200),new.sim$data$censored.time)
  expect_equal(new.sim$data$Id,sim.dropout$data$Id)
  expect_equal(new.sim$data$arm,sim.dropout$data$arm)
  expect_equal(new.sim$data$actual.events,sim.dropout$data$actual.events)
  
  expect_true(all(new.sim$data$observed.events >= sim.dropout$data$observed.events))
  expect_equal(new.sim$data$observed.events, vapply(new.sim$event.times,length,FUN.VALUE = numeric(1)))
  
  invisible(mapply(function(x,y,z){expect_equal(c(x,y),z)},
         x=sim.dropout$event.times,
         y=impute$imputed.values[,4]$newevent.times,
         z=new.sim$event.times))
  
  #Note no dropouts as far as new.sim is concerned
  expect_equal(c(0,0),summary(new.sim)$number.dropouts)
  
})

test_that("ImputeSim.Simfit",{
  
})
