context("ImportTests")

test_that("simple_invalid_args",{
  
  covar.df <- data.frame(Id=1:6,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <- list(c(25.6,100,121,200,225),c(100,110),c(55),numeric(0),150,45)
  
  #check valid
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="complete",study.time=365),
               regexp=NA)
  
  
  #dejaData invalid class
  expect_error(ImportSim(dejaData=covar.df,event.times=event.times,
                         status="complete",study.time=365))
  
  #invalid study.time
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="complete",study.time=c(10,20)))
  
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="complete",study.time=-4))
  
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="complete",study.time=0))
  
  #invalid status
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="meh",study.time=365))
})



test_that("event.times",{
  
  covar.df <- data.frame(Id=1:6,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <-
      list(c(25.6,100,121,200,225),
           c(100,110),
           c(55),
           numeric(0),
           150,45)
  
  #wrong number of event.times
  e <- event.times
  e[[7]] <- numeric(0)
  
  expect_error(ImportSim(dejaData=dejaData,
                         event.times=e,
                         status="complete",
                         study.time=365))
  
  #event.times not sorted
  e <- event.times
  e[[6]] <- c(60,20)
  expect_error(ImportSim(dejaData=dejaData,event.times=e,
                         status="complete",study.time=365))
  
  #non-numeric argument
  e <- event.times
  e[[1]] <- logical(0)
  expect_error(ImportSim(dejaData=dejaData,event.times=e,
                         status="complete",study.time=365))
  
  #negative event.time
  e <- event.times
  e[[5]] <- c(-10,10)
  expect_error(ImportSim(dejaData=dejaData,event.times=e,
                         status="complete",study.time=365))
  
  #event time > study.time
  e <- event.times
  e[[2]] <- c(5,10,100,900)
  expect_error(ImportSim(dejaData=dejaData,event.times=e,
                         status="complete",study.time=365))

  ## event time > study time and allow.beyond.study
  expect_error(ImportSim(dejaData=dejaData,
                         event.times=e,
                         status="complete",
                         study.time=365,
                         allow.beyond.study=TRUE),
               regexp=NA)
})

test_that("Invalid_args_for_complete",{
  
  covar.df <- data.frame(Id=1:6,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <- list(
      c(25.6,100,121,200,225),
      c(100,110),
      c(55),
      numeric(0),
      150,
      45)
  
  #cannot use censored.time
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="complete",study.time=365,censored.time=rep(365,6)))
  
  #can use actual.events
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="complete",study.time=365,actual.events=rep(20,6)))
  
})

test_that("valid_args_complete",{
  covar.df <- data.frame(Id=10:15,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0),Z2=c(4,5,1,3,2,4))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <- list(c(25.6,100,121,200,225),c(100,110),c(55),numeric(0),150,45)
  
  sim <- ImportSim(dejaData=dejaData,event.times=event.times,status="complete",study.time=365)
  
  expect_equal("SingleSim",class(sim))
  expect_equal("complete",sim$status)
  expect_equal(event.times,sim$event.times)
  expect_true(is.null(sim$subject.rates))
  expect_true(is.null(sim$dropout.mechanism))
  expect_true(is.null(sim$impute.mechanism))
  expect_equal(365,sim$study.time)
  expect_true(is.null(sim$event.rates))
  expect_true(is.null(sim$dispersions))
  
  expect_equal(covar.df$Z,sim$data$Z)
  expect_equal(covar.df$Z2,sim$data$Z2)
  expect_equal(sim$data$observed.events,sim$data$actual.events)
  expect_equal(c(5,2,1,0,1,1),sim$data$observed.events)
  expect_equal(rep(365,6),sim$data$censored.time)
  expect_equal(10:15,sim$data$Id)
})

test_that("invalid_censored_time",{
  covar.df <- data.frame(Id=10:15,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0),Z2=c(4,5,1,3,2,4))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <- list(c(25.6,100,121,200,225),c(100,110),c(55),numeric(0),150,45)
  
  #event after censoring
  censored.time <- c(200,365,365,365,365,365)
  
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout",study.time=365,censored.time=censored.time))
  
  
  #No censored.time argument
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout",study.time=365))
  
  
  censored.time <- c(365,365,365,-5,365,365)
  #negative time
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout",study.time=365,censored.time=censored.time))
  
  #censored > study.time
  censored.time <- rep(400,6)
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout",study.time=365,
                         censored.time=censored.time))

  ## allow censored > study.time if allow.beyond.study is TRUE
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout", study.time=365,
                         censored.time=censored.time,
                         allow.beyond.study=TRUE),
               regexp=NA)
  
  #wrong length
  censored.time <- rep(300,7)
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout",study.time=365,censored.time=censored.time))
  
  #not number
  censored.time <- rep(TRUE,6)
  expect_error(ImportSim(dejaData=dejaData,event.times=event.times,
                         status="dropout",study.time=365,censored.time=censored.time))
})


test_that("invalid_actual.events",{
  covar.df <- data.frame(Id=10:15,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0),Z2=c(4,5,1,3,2,4))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <- list(c(25.6,100,121,200,225),c(100,110),c(55),numeric(0),150,45)
  
  censored.time <- c(365,250,365,30,200,50)

  #wrong length
  actual.events <- rep(10,5)
  expect_error(ImportSim(dejaData=dejaData,
                         event.times=event.times,
                         status="dropout",
                         study.time=365,
                         censored.time=censored.time,
                         actual.events=actual.events))
  
  #not integers
  actual.events <- rep(10.5,6)
  expect_error(ImportSim(dejaData=dejaData,
                         event.times=event.times,
                         status="dropout",
                         study.time=365,
                         censored.time=censored.time,
                         actual.events=actual.events))
  
  
  #< observed.events
  actual.events <- c(4,3,1,1,2,5)
  expect_error(ImportSim(dejaData=dejaData,
                         event.times=event.times,
                         status="dropout",
                         study.time=365,
                         censored.time=censored.time,
                         actual.events=actual.events))
  
})


test_that("valid_dropout",{
  #first no actual events
  covar.df <- data.frame(Id=10:15,arm=c(rep(0,3),rep(1,3)),Z=c(0,1,1,0,1,0),Z2=c(4,5,1,3,2,4))
  dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id")
  event.times <- list(c(25.6,100,121,200,225),c(100,110),c(55),numeric(0),150,45)
  
  censored.time <- c(365,250,365,30,200,50)
  
  sim <- ImportSim(dejaData=dejaData,event.times=event.times,
                   status="dropout",study.time=365,censored.time=censored.time)
  
  #only tests specifically to do with dropout are here, see complete case
  #for more general tests
  expect_equal("dropout",sim$status)
  expect_equal(censored.time,sim$data$censored.time)
  expect_equal(rep(as.numeric(NA),6),sim$data$actual.events)
  expect_equal("User imported dropout data",sim$dropout.mechanism)
  
  #now use actual.events
  sim <- ImportSim(dejaData=dejaData,event.times=event.times,
                   status="dropout",study.time=365,censored.time=censored.time,
                   actual.events=c(5,NA,1,3,1,NA))
  expect_equal(c(5,NA,1,3,1,NA),sim$data$actual.events)
  expect_equal(c(5,2,1,0,1,1),sim$data$observed.events)
})
