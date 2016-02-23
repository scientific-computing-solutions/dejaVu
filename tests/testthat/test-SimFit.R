context("SingleSimFit")

test_that("GetDefaultFormula",{
  
  expect_error(GetDefaultFormula(c(TRUE,TRUE)))
  expect_error(GetDefaultFormula("x"))
  
  expect_equal(as.character(GetDefaultFormula(TRUE)),as.character(formula(observed.events~arm+offset(log(censored.time)))))
  expect_equal(as.character(GetDefaultFormula(FALSE)),as.character(formula(observed.events~offset(log(censored.time)))))
  
})

test_that("Simfit.SingleSim",{
  #not testing formula argument
  
  set.seed(143)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  sim.dropout <- SimulateDropout(sim, 
                    drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  
  expect_error(Simfit("we"))
  expect_error(Simfit(sim.dropout,family="gaussian"))
  
  expect_error(Simfit(sim.dropout,family="poisson",equal.dispersion=FALSE))
})


test_that("Simfit_diff_dispersions",{
  set.seed(143)

  sim <- SimulateComplete(study.time=365, 
                        number.subjects=50, 
                        event.rates=c(0.01,0.005),
                        dispersions=0.25) 

  sim.dropout <- SimulateDropout(sim, 
                               drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  fit.diff.disp <- Simfit(sim.dropout,equal.dispersion=FALSE)
  expect_error(summary(fit.diff.disp))
  
  expect_equal("SingleSimFit",class(fit.diff.disp))
  
  expect_equal("SingleSim",class(fit.diff.disp$singleSim))
  expect_equal(fit.diff.disp$singleSim$data,sim.dropout$data)
  
  expect_false(fit.diff.disp$impute.parameters$equal.dispersion)
  
  mymod.1 <- glm.nb(observed.events ~ offset(log(censored.time)) ,data=sim.dropout$data[sim.dropout$data$arm==0,])
  mymod.2 <- glm.nb(observed.events ~ offset(log(censored.time)) ,data=sim.dropout$data[sim.dropout$data$arm==1,])
  #hack calls so they match
  mymod.1$call <- fit.diff.disp$model[[1]]$call
  mymod.2$call <- fit.diff.disp$model[[2]]$call
  expect_equal(mymod.1,fit.diff.disp$model[[1]])
  expect_equal(mymod.2,fit.diff.disp$model[[2]])
  
  expect_equal(c(mymod.1$theta,mymod.2$theta), fit.diff.disp$impute.parameters$gamma)
  mu <- matrix(rep(exp(c(coefficients(mymod.1),coefficients(mymod.2))),100),ncol=2,byrow=TRUE)
  names(mu) <- NULL
  expect_equal(mu, fit.diff.disp$impute.parameters$mu)
})


test_that("sim_fit_equal_dispersion",{
  set.seed(3143)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  sim.dropout <- SimulateDropout(sim, 
                                 drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  fit <- Simfit(sim.dropout) #equal.disperison is default
    
  expect_equal(100,numberSubjects(fit))
  expect_equal(fit$singleSim$data,sim.dropout$data)
  
  mymod <- glm.nb(observed.events ~arm + offset(log(censored.time)),data=sim.dropout$data)
  #hack call to make them equal
  mymod$call <- fit$model$call
  expect_equal(mymod,fit$model)
  
  expect_true(fit$impute.parameters$equal.dispersion)
  expect_equal(rep(mymod$theta,2),fit$impute.parameters$gamma )
  
  #checking mu
  mu <- mu <- matrix(rep(exp(mymod$coefficients[1])*c(1,exp(mymod$coefficients[2])),100),ncol=2,byrow=TRUE)
  names(mu) <- NULL
  expect_equal(mu,fit$impute.parameters$mu )
  
})


test_that("simfit_poisson_and_qpoisson",{
  set.seed(1435)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  fit.pois <- Simfit(sim,family="poisson")
  
  params <- list(equal.dispersion=TRUE)
  class(params) <- "ImputeParameters"
  
  expect_equal(params,fit.pois$impute.parameters)
  mod <- glm(observed.events ~ arm + offset(log(censored.time)), family="poisson",data=sim$data)
  #hack call so they match
  mod$call <- fit.pois$model$call
  expect_equal(mod,fit.pois$model)
  
  fit.qpois <- Simfit(sim,family="quasipoisson")
  
  expect_equal(params,fit.qpois$impute.parameters)
  mod <- glm(observed.events ~ arm + offset(log(censored.time)), family="quasipoisson",data=sim$data)
  #hack call so they match
  mod$call <- fit.qpois$model$call
  expect_equal(mod,fit.qpois$model)
  
})

test_that("summary.SingleSimFit",{
  set.seed(3143)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  sim.dropout <- SimulateDropout(sim, 
                                 drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  
  fit <- Simfit(sim.dropout)
  
  expect_error(summary(fit,CI.limit=-1))
  expect_error(summary(fit,CI.limit=1))
  
  s <- summary(fit)
  
  expect_equal("summary.SingleSimFit",class(s))
  expect_equal(summary(fit$model), s$model.summary)
  expect_equal(0.95,s$CI.limit)
  expect_equal(summary(fit$singleSim)$number.dropouts, s$dropout)
  tr <- coefficients(fit$model)
  names(tr) <- NULL
  expect_equal(exp(tr[2]),s$treatment.effect)
  expect_equal(fit$impute.parameters$mu[1,], s$rate.estimate)
  expect_equal(fit$model$df.residual,s$df)
  expect_equal("dropout",s$datastatus)
  expect_equal(1/fit$impute.parameters$gamma[1], s$dispersion)
  expect_equal(coefficients(summary(fit$model))[2,2],s$se)
  expect_equal(2*pnorm(log(s$treatment.effect)/s$se),s$pval)
  
  s2 <- summary(fit,CI.limit = 0.75)
  expect_true(s2$CI[1] > s$CI[1])
  expect_true(s2$CI[2] < s$CI[2])
  
})

test_that("summary.SingleSimFit_for_poisson",{
  set.seed(33143)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  sim.dropout <- SimulateDropout(sim, 
                                 drop.mechanism=ConstantRateDrop(rate=0.0025,var=1))
  
  fit <- Simfit(sim.dropout,family="poisson")
  
  s <- summary(fit)
  
  expect_equal("summary.SingleSimFit",class(s))
  expect_equal(summary(fit$model), s$model.summary)
  expect_equal(0.95,s$CI.limit)
  expect_equal(summary(fit$singleSim)$number.dropouts, s$dropout)
  tr <- coefficients(fit$model)
  names(tr) <- NULL
  expect_equal(exp(tr[2]),s$treatment.effect)
  expect_equal(exp(tr[1])*c(1,exp(tr[2])), s$rate.estimate)
  expect_equal(fit$model$df.residual,s$df)
  expect_equal("dropout",s$datastatus)
  expect_equal(numeric(0), s$dispersion)
  expect_equal(coefficients(summary(fit$model))[2,2],s$se)
  expect_equal(2*pnorm(log(s$treatment.effect)/s$se),s$pval)
})

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
  
  expect_equal(fit$impute.parameters,imputed$impute.parameters)
  expect_equal(impute.mechanism,imputed$impute.mechanism)
  expect_equal(sim.dropout$data,imputed$singleSim$data)
  
  expect_equal(1,ncol(imputed$imputed.values))
  expect_equal(2,nrow(imputed$imputed.values))
  
  expect_equal(rep(365,200),imputed$imputed.values[,1]$new.censored.times)
  
  retVal <- ifelse(sim.dropout$data$arm==1 | sim.dropout$data$censored.time==365,NA,365)
  retVal <- lapply(retVal,function(x){if(is.na(x)) numeric(0) else x})
  
  expect_equal(retVal,imputed$imputed.values[,1]$newevent.times)
  
})

