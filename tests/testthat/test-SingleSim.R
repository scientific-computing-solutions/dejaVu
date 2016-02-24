context("SingleSim")

test_that("getDispersions",{
  expect_equal(c(0.5,0.5,0.5,0.5),getDispersions(c(0,0,0,0),c(0.5,0.3)))
  expect_equal(c(0.5,0.3,0.3,0.5),getDispersions(c(0,1,1,0),c(0.5,0.3)))
  expect_equal(c(0.3,0.3,0.3,0.3),getDispersions(c(0,1,1,0),0.3))
})

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
  
  expect_error(SimulateComplete(study.time=10,dejaData = "hello",dispersions=5))
  dejaData <- MakeDejaData(data=data.frame(Id=1:10,arm=c(0,rep(1,9))),Id="Id",arm="arm")
  
  #No rate in dejaData
  expect_error(SimulateComplete(study.time=10,dejaData =dejaData,dispersions=5))
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


test_that("numberSubjects_subjectsperarm",{
  sim <- SimulateComplete(study.time=12,number.subjects=c(2,4),event.rates=0.1,dispersions=c(0,0.5))
  expect_equal(c(2,4),subjectsPerArm(sim))
  expect_equal(6,numberSubjects(sim))
})

test_that("SimComplete_with_deja_matches_without",{
  set.seed(20)
  sim <- SimulateComplete(study.time=10,number.subjects=2,event.rates=c(0.5,0.6),dispersions=0.5)
  
  set.seed(20)
  dejaData <- MakeDejaData(data=data.frame(Id=1:4,arm=c(0,0,1,1),rate=c(0.5,0.5,0.6,0.6)),
                           Id="Id",arm="arm",rate="rate")
  simdeja <- SimulateComplete(study.time=10,dejaData=dejaData,dispersions=0.5)
  
  #the event.rates are different as they are not included in simdeja so remove them before
  #checking for equality
  sim$event.rates <- NULL
  simdeja$event.rates <- NULL
  expect_equal(sim,simdeja)
  
})


test_that("Simcomplete_with_deja",{
  set.seed(20)
  dejaData <- MakeDejaData(data=data.frame(Z=c(0,1,1,0),id=1:4,Arm=c(0,0,1,1),rate=c(0.5,0,0.5,0)),
                           Id="id",arm="Arm",rate="rate")
  simdeja <- SimulateComplete(study.time=10,dejaData=dejaData,dispersions=0)
  
  expect_equal(c("Id","arm","censored.time","observed.events","actual.events","Z"),colnames(simdeja$data))
  expect_equal(c(0,1,1,0),simdeja$data$Z)
  
  #rate zero gives us zero
  expect_equal(c(0,0),simdeja$data$observed.events[c(2,4)])
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


test_that("SingleSimSummary",{
  set.seed(400)
  sim <- SimulateComplete(study.time=12,number.subjects=c(2,4),event.rates=0.1,dispersions=c(0,0.5))
  
  #for complete dataset
  s <- summary(sim)
  expect_equal("complete",s$status)
  expect_equal(12,s$study.time)
  expect_equal(c(2,4),s$number.subjects)
  expect_equal(sum(sim$data[sim$data$arm==0,]$observed.events), s$total.events[1])
  expect_equal(sum(sim$data[sim$data$arm==1,]$observed.events), s$total.events[2])
  expect_equal(c(0,0),s$number.dropouts)
  expect_equal(sum(sim$data[sim$data$arm==0,]$censored.time), s$time.at.risk[1])
  expect_equal(sum(sim$data[sim$data$arm==1,]$censored.time), s$time.at.risk[2])
  expect_equal(s$empirical.rates,s$total.events/s$time.at.risk)
  
  dummy.dropout <- CreateNewDropoutMechanism(type="MAR",text="hello",
                                             GetDropTime=function(event.times,data){
                                               return(data$censored.time)
                                             })
  
  dropout.sim <- SimulateDropout(sim,dummy.dropout)
  #hack dropout sim to get c(1,0) dropouts
  dropout.sim$data[1,]$censored.time <- 11
  dropout.sim$data[1,]$observed.events <- 1
  
  s <- summary(dropout.sim)
  expect_equal(c(1,0),s$number.dropouts)
  expect_equal(sum(dropout.sim$data[dropout.sim$data$arm==0,]$observed.events), s$total.events[1])
  expect_equal(sum(dropout.sim$data[dropout.sim$data$arm==1,]$observed.events), s$total.events[2])
  expect_equal(sum(dropout.sim$data[dropout.sim$data$arm==0,]$censored.time), s$time.at.risk[1])
  expect_equal(sum(dropout.sim$data[dropout.sim$data$arm==1,]$censored.time), s$time.at.risk[2])
  
})
