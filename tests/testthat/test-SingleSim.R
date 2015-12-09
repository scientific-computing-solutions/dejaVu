context("SingleSim")

test_that("Simulate_Complete_Invalid_Args",{
  expect_error(SimulateComplete(study.time=-1,number.subjects=50,event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=0,number.subjects=50,event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=c(5,6),number.subjects=50,event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=10,number.subjects="x2",event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=10,number.subjects=c(2,5,10),event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=10,number.subjects=c(0,5),event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=10,number.subjects=c(1.4,5),event.rates=0.5,dispersions=0.5))
  expect_error(SimulateComplete(study.time=10,number.subjects=2,event.rates=c(-0.5,5),dispersions=0.5))
  expect_error(SimulateComplete(study.time=10,number.subjects=2,event.rates=0.5,dispersions=-5))
})

test_that("CompleteSim_creation_as_expected",{
  #looking at non stochastic parts
  set.seed(20)
  sim <- SimulateComplete(study.time=10,number.subjects=2,event.rates=0.5,dispersions=0.5)
  
  expect_equal("SingleSim",class(sim))
  expect_equal("complete",sim$status)
  expect_equal(4,numberSubjects(sim))
  expect_equal(4,nrow(sim$data))
  expect_true(is.null(sim$dropout.mechanism))
  expect_true(is.null(sim$impute.mechanism))
  expect_equal(10,sim$study.time)
  expect_equal(0.5,sim$dispersions)
  expect_equal(4,length(sim$subject.rates))
  expect_equal(4,length(sim$event.times))
  expect_equal(10,sim$study.time)
  data <- sim$data
  expect_equal(data$observed.events, vapply(sim$event.times,length,FUN.VALUE = numeric(1)))
  expect_equal(data$observed.events,data$actual.events)
  expect_equal(as.factor(c(0,0,1,1)),data$arm)
  expect_equal(1:4,data$Id)
  expect_equal(rep(10,4),data$censored.time)
  
  #dispersion can equal 0
  sim <- SimulateComplete(study.time=12,number.subjects=c(2,4),event.rates=0.1,dispersions=c(0,0.5))
  expect_equal("SingleSim",class(sim))
  expect_equal("complete",sim$status)
  expect_equal(6,numberSubjects(sim))
  expect_equal(6,nrow(sim$data))
  expect_true(is.null(sim$dropout.mechanism))
  expect_true(is.null(sim$impute.mechanism))
  expect_equal(12,sim$study.time)
  expect_equal(c(0,0.5),sim$dispersions)
  expect_equal(6,length(sim$subject.rates))
  expect_equal(6,length(sim$event.times))
  expect_equal(12,sim$study.time)
  data <- sim$data
  expect_equal(data$observed.events, vapply(sim$event.times,length,FUN.VALUE = numeric(1)))
  expect_equal(data$observed.events,data$actual.events)
  expect_equal(as.factor(c(0,0,1,1,1,1)),data$arm)
  expect_equal(1:6,data$Id)
  expect_equal(rep(12,6),data$censored.time)
  expect_true(all(unlist(c(sim$event.times))<=12))
  invisible(lapply(1:6,function(i)expect_true(all(sim$event.times[[i]]==sort(sim$event.times[[i]])))))
})

test_that("SimulateDropout",{
  #non stochastic parts behave as expected
  set.seed(9)
  sim <- SimulateComplete(study.time=12,number.subjects=c(2,4),event.rates=0.1,dispersions=c(0,0.5))
  dummy.dropout <- CreateNewDropoutMechanism(type="MAR",text="hello",
                                             GetDropTime=function(event.times,data){
                                               return(data$censored.time/2)
                                             })
  
  expect_error(SimulateDropout(3,dummy.dropout))
  expect_error(SimulateDropout(sim,data.frame(x=c(1,2,3))))
  
  sim2 <- sim
  sim2$status <- "imputed"
  expect_error(SimulateDropout(sim2,dummy.dropout))
  
  sim.dropout <- SimulateDropout(sim,dummy.dropout)
  expect_equal("SingleSim",class(sim.dropout))
  expect_equal("dropout",sim.dropout$status)
  
  expect_equal(sim.dropout$data$actual.events,sim$data$observed.events)
  expect_equal(rep(6,6),sim.dropout$data$censored.time)
  
  expect_equal(sim.dropout$data$observed.events, vapply(sim.dropout$event.times,length,FUN.VALUE = numeric(1)))
  
  expect_true(all(unlist(c(sim.dropout$event.times))<=6))
  
  invisible(lapply(1:6,function(i){ et <- sim$event.times[[i]]
                                    et <- et[et<=6]
                                    expect_equal(et,sim.dropout$event.times[[i]])}))
  
  expect_equal(sim.dropout$dropout.mechanism$text,dummy.dropout$text)
})

