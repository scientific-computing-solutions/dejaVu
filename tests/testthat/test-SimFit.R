context("SingleSimFit")

test_that("GetModelFormula",{
  
  expect_error(GetModelFormula(c(TRUE,TRUE),covar=NULL))
  expect_error(GetModelFormula("x",covar=NULL))
  
  expect_equal(as.character(GetModelFormula(TRUE,covar=NULL)),as.character(formula(observed.events~arm+offset(log(censored.time)))))
  expect_equal(as.character(GetModelFormula(FALSE,covar=NULL)),as.character(formula(observed.events~offset(log(censored.time)))))
  
})

test_that("Simfit.SingleSim",{
  #not testing formula argument here
  
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
  
  expect_false(fit.diff.disp$equal.dispersion)
  
  mymod.1 <- glm.nb(observed.events ~ offset(log(censored.time)) ,data=sim.dropout$data[sim.dropout$data$arm==0,])
  mymod.2 <- glm.nb(observed.events ~ offset(log(censored.time)) ,data=sim.dropout$data[sim.dropout$data$arm==1,])
  #hack calls so they match
  mymod.1$call <- fit.diff.disp$model[[1]]$call
  mymod.2$call <- fit.diff.disp$model[[2]]$call
  expect_equal(mymod.1,fit.diff.disp$model[[1]])
  expect_equal(mymod.2,fit.diff.disp$model[[2]])
  
  expect_warning(expect_equal(c(mymod.1$theta,mymod.2$theta), fit.diff.disp$genCoeff.function(FALSE)$gamma))
  mu <- matrix(rep(exp(c(coefficients(mymod.1),coefficients(mymod.2))),100),ncol=2,byrow=TRUE)
  names(mu) <- NULL
  expect_warning(expect_equal(mu, fit.diff.disp$genCoeff.function(FALSE)$mu))
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
  
  expect_true(fit$equal.dispersion)
  expect_warning(expect_equal(rep(mymod$theta,2),fit$genCoeff.function(FALSE)$gamma))
  
  #checking mu
  mu <- matrix(rep(exp(mymod$coefficients[1])*c(1,exp(mymod$coefficients[2])),100),ncol=2,byrow=TRUE)
  names(mu) <- NULL
  expect_warning(expect_equal(mu,fit$genCoeff.function(FALSE)$mu))
  
})


test_that("simfit_poisson_and_qpoisson",{
  set.seed(1435)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  fit.pois <- Simfit(sim,family="poisson")
  
  expect_true(fit.pois$equal.dispersion)
  expect_true(is.null(fit.pois$genCoeff.function))
  mod <- glm(observed.events ~ arm + offset(log(censored.time)), family="poisson",data=sim$data)
  #hack call so they match
  mod$call <- fit.pois$model$call
  expect_equal(mod,fit.pois$model)
  
  fit.qpois <- Simfit(sim,family="quasipoisson")
  
  expect_true
  expect_true(is.null(fit.qpois$genCoeff.function))
  mod <- glm(observed.events ~ arm + offset(log(censored.time)), family="quasipoisson",data=sim$data)
  #hack call so they match
  mod$call <- fit.qpois$model$call
  expect_equal(mod,fit.qpois$model)
  
})


test_that("SimFit_covar",{
  
  set.seed(20)
  dejaData <- MakeDejaData(data=data.frame(Z=rep(c(0,1),20),
                                           id=1:40,
                                           Arm=c(rep(0,20),rep(1,20)),
                                           rate=c(rep(0.5,40))),
                           Id="id",arm="Arm",rate="rate")
  simdeja <- SimulateComplete(study.time=10,dejaData=dejaData,dispersions=0)
  
  expect_error(Simfit(simdeja,covar=~arm))
  expect_error(Simfit(simdeja,covar="hello"))
  expect_error(Simfit(simdeja,covar=~arm*Z))
  
  sf <- Simfit(simdeja,covar=~Z)
  mod <- glm.nb(observed.events ~ arm + Z + offset(log(censored.time)),data=simdeja$data)
  expect_equal(sf$model$coefficients,mod$coefficients)
  expect_equal(vcov(sf$model),vcov(mod))
  #apart from call they should be the same
  mod$call <- sf$model$call
  expect_equal(sf$model,mod)
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
  expect_warning(expect_equal(fit$genCoeff.function(FALSE)$mu[1,], s$rate.estimate))
  expect_equal(fit$model$df.residual,s$df)
  expect_equal("dropout",s$datastatus)
  expect_warning(expect_equal(1/fit$genCoeff.function(FALSE)$gamma[1], s$dispersion))
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

