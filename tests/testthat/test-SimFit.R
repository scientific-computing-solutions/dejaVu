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
  p <- exp(c(coefficients(mymod.1),coefficients(mymod.2)))
  names(p) <- NULL
  expect_equal(p, fit.diff.disp$impute.parameters$p)
})

test_that("simfit_poisson_and_qpoisson",{
  set.seed(1435)
  
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=50, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  fit.pois <- Simfit(sim,family="poisson")
  
  expect_equal(list(equal.dispersion=TRUE),fit.pois$impute.parameters)
  mod <- glm(observed.events ~ arm + offset(log(censored.time)), family="poisson",data=sim$data)
  #hack call so they match
  mod$call <- fit.pois$model$call
  expect_equal(mod,fit.pois$model)
  
  fit.qpois <- Simfit(sim,family="quasipoisson")
  
  expect_equal(list(equal.dispersion=TRUE),fit.qpois$impute.parameters)
  mod <- glm(observed.events ~ arm + offset(log(censored.time)), family="quasipoisson",data=sim$data)
  #hack call so they match
  mod$call <- fit.qpois$model$call
  expect_equal(mod,fit.qpois$model)
  
})