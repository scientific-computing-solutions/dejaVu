context("Impute")

test_that("Impute_valid_args",{
  set.seed(33143)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=100, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  sim.dropout <- SimulateDropout(sim, 
                                 drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  
  fit <- Simfit(sim.dropout)
  
  impute.mechanism <- weighted_j2r(0.5)
  
  expect_error(Impute(fit="sd",impute.mechanism=impute.mechanism,N=10))
  expect_error(Impute(fit=fit,impute.mechanism=c(1,2,3),N=10))
  expect_error(Impute(fit=fit,impute.mechanism=impute.mechanism,N=0))
  expect_error(Impute(fit=fit,impute.mechanism=impute.mechanism,N=11.5))
  expect_error(Impute(fit=fit,impute.mechanism=impute.mechanism,N=c(1,0.5)))
})

test_that("Impute_general_mechanism",{
  
  set.seed(33143)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=100, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  sim.dropout <- SimulateDropout(sim, 
                                 drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  
  fit <- Simfit(sim.dropout)
  impute.mechanism <- CreateNewImputeMechanism(name="my.imp",cols.needed=c("censored.time","arm"),parameters=list(W=1),
                                               impute=function(fit){
                                                 df <- fit$singleSim$data
                                                 number.of.subjects <- numberSubjects(fit$singleSim)
                                                 study.time <- fit$singleSim$study.time
                                                 
                                                 newevent.times <- lapply(1:number.of.subjects,function(i){
                                                   time.left <- fit$singleSim$study.time - df$censored.time[i]
                                                   if(df$arm[i]==1|| time.left==0) return(numeric(0))
                                                   new.event.times <- fit$singleSim$study.time  
                                                   
                                                 })
                                                 
                                                 return(list(new.censored.times=rep(fit$singleSim$study.time,numberSubjects((fit))),
                                                             newevent.times=newevent.times))
                                                 
                                               })
  
  imputed <- Impute(fit,impute.mechanism,N=1)
  
  expect_equal("ImputeSim",class(imputed))
  expect_equal(summary(sim.dropout)$number.dropouts,imputed$dropout)
  
  expect_equal(impute.mechanism,imputed$impute.mechanism)
  expect_equal(sim.dropout$data,imputed$singleSim$data)
  
  expect_equal(1,ncol(imputed$imputed.values))
  expect_equal(2,nrow(imputed$imputed.values))
  
  expect_equal(rep(365,200),imputed$imputed.values[,1]$new.censored.times)
  
  retVal <- ifelse(sim.dropout$data$arm==1 | sim.dropout$data$censored.time==365,NA,365)
  retVal <- lapply(retVal,function(x){if(is.na(x)) numeric(0) else x})
  
  expect_equal(retVal,imputed$imputed.values[,1]$newevent.times)
  
})

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
  set.seed(202)
  
  sim <- SimulateComplete(study.time=365,number.subjects=100,
                          event.rates=c(0.01,0.005),dispersions=0.25)
  
  sim.dropout <- SimulateDropout(sim,drop.mechanism=ConstantRateDrop(rate=0.0025)) 
  fit <- Simfit(sim.dropout,equal.dispersion=TRUE)
  
  impute <- Impute(fit,impute.mechanism = weighted_j2r(trt.weight=0),N=8)
  
  expect_error(Simfit(impute,equal.dispersion = FALSE))
  
  checkfits <- function(family){
    fits <- Simfit(impute,family=family)
    expect_equal("ImputeSimFit",class(fits))
    expect_equal(fits$imputeSim,impute)
    invisible(lapply(fits$summaries,function(x){expect_equal("summary.SingleSimFit",class(x))}))
    expect_equal(fits$summaries[[5]],summary(Simfit(GetImputedDataSet(impute,5),family=family)))
    expect_equal(8,length(fits$summaries))
  }
    
  checkfits(family="negbin")
  checkfits(family="poisson")
  checkfits(family="quasipoisson")
  
  
})

test_that("as.data.frame.ImputeSimfit",{
  
  set.seed(1202)
  
  sim <- SimulateComplete(study.time=365,number.subjects=100,
                          event.rates=c(0.07,0.05),dispersions=c(0,0.25))
  
  sim.dropout <- SimulateDropout(sim,drop.mechanism=ConstantRateDrop(rate=0.0025)) 
  fit <- Simfit(sim.dropout,equal.dispersion=TRUE)
  
  impute <- Impute(fit,impute.mechanism = weighted_j2r(trt.weight=0),N=18)
  fits <- Simfit(impute)
  
  my.df <- as.data.frame(fits)
  expected.cols <- c("imputeID","control.rate","active.rate","treatment.effect","se","pval","dispersion")
  
  expect_equal(expected.cols,colnames(my.df))
  expect_equal(18,nrow(my.df))
  
  expect_equal(1:18,my.df$imputeID)
  
  expect_equal(my.df$se[12],fits$summaries[[12]]$se)
  expect_equal(my.df$pval[15],fits$summaries[[15]]$pval)
  expect_equal(my.df$dispersion[5],fits$summaries[[5]]$dispersion)
  expect_equal(my.df$control.rate[7],fits$summaries[[7]]$rate.estimate[1])
  expect_equal(my.df$active.rate[7],fits$summaries[[7]]$rate.estimate[2])
  
  expect_error(Simfit(impute,family=poisson))
  
  fits <- Simfit(impute,family="poisson")
  my.df <- as.data.frame(fits)
  expect_true(all(is.na(my.df$dispersion)))
  
})


test_that("summary.ImputeSimfit.fails.if.1.dataset",{
  set.seed(1202)
  
  sim <- SimulateComplete(study.time=365,number.subjects=100,
                          event.rates=c(0.07,0.05),dispersions=c(0,0.25))
  
  sim.dropout <- SimulateDropout(sim,drop.mechanism=ConstantRateDrop(rate=0.0025)) 
  fit <- Simfit(sim.dropout,equal.dispersion=TRUE)
  
  impute <- Impute(fit,impute.mechanism = weighted_j2r(trt.weight=0),N=1)
  fits <- Simfit(impute)
  
  expect_error(summary(fits))
  
})


test_that("summary.ImputeSimfit",{
  set.seed(15202)
  
  sim <- SimulateComplete(study.time=365,number.subjects=100,
                          event.rates=c(0.07,0.05),dispersions=c(0,0.25))
  
  sim.dropout <- SimulateDropout(sim,drop.mechanism=ConstantRateDrop(rate=0.0025)) 
  fit <- Simfit(sim.dropout,equal.dispersion=TRUE)
  
  impute <- Impute(fit,impute.mechanism = weighted_j2r(trt.weight=0),N=18)
  fits <- Simfit(impute)
  
  s <- summary(fits)
  expect_equal("summary.ImputeSimFit",class(s))
  expect_equal(s$dropout,fits$imputeSim$dropout)
  dispersions <- vapply(fits$summaries,function(x)x$dispersion,numeric(1))
  
  expect_equal(mean(dispersions),s$dispersion)
  
  s <- summary(Simfit(impute,family="quasipoisson"))
  expect_true(all(is.na(s$dispersion)))
})

test_that("Rubins.formula",{
  treatment.effect <- c(3,4,5,6)
  ses <- c(1,2,3,6)
  original.df <- 40
  N <- 4
  
  retVal <- .rubinsformula(treatment.effect,ses,original.df,N)
  
  expect_equal(exp(mean(log(treatment.effect))),retVal$treatment.effect)
  se <- sqrt((5/4)*var(c(log(3:6))) + mean(ses^2))
  expect_equal(se,retVal$se)
  
  df <- 3*(1+mean(ses^2)/(5*var(c(log(3:6)))/4))^2
  expect_equal(df,retVal$df)
  
  v <- (40*41/43)*(1-(5/4)*(var(c(log(3:6)))/se^2))
  adjusted.df <- 1/(1/df + 1/v )
  expect_equal(adjusted.df,retVal$adjusted.df)
})
