context("ImputeMechanisms")

#Note no validation of GetDropTime function

test_that("CreateNewImputeMechanism_invalid_args",{
  expect_error(CreateNewImputeMechanism(name=c(1,2,3),impute=function(){})) 
  expect_error(CreateNewImputeMechanism(name="hello",impute=45)) 
  expect_error(CreateNewImputeMechanism(name="hello",impute=function(){},parameters = c(3,4)))
})

test_that("DummyImputeMechanism",{
  
  imp <- CreateNewImputeMechanism(name="my.imp",cols.needed="censored.time",parameters=list(W=1),
                                  impute=function(fit){
                                    return(list(new.censored.times=rep(fit$singleSim$study.time,numberSubjects((fit))),
                                                newevent.times=rep(numeric(0),numberSubjects(fit))))
                                    
                                  })
  
  expect_equal("ImputeMechanism",class(imp))
  expect_equal("my.imp",imp$name)
  expect_equal("censored.time",imp$cols.needed)
  expect_equal(list(W=1),imp$parameters)

})

test_that("weightj2r",{
  expect_error(weighted_j2r(trt.weight=-1))
  expect_error(weighted_j2r(trt.weight=1.5))
  expect_error(weighted_j2r(trt.weight="dfge"))
  expect_error(weighted_j2r(trt.weight=0.5,delta=c(4,5)))
  expect_error(weighted_j2r(trt.weight=0,delta=c(4,5)))
  expect_error(weighted_j2r(trt.weight=1,delta=c(4,5,6)))
  expect_error(weighted_j2r(trt.weight=1,delta=c(-5,5)))
  expect_error(weighted_j2r(trt.weight=1,delta=c(1.5,"sdf")))
  expect_error(weighted_j2r(trt.weight=1,delta=c(0,1)))
  
  x <- weighted_j2r(trt.weight=0.5)
  expect_equal("ImputeMechanism",class(x))
  expect_equal("weighted_j2r",x$name)
  expect_equal(c("censored.time","observed.events","arm"),x$cols.needed)
  expect_equal(list(trt.weight=0.5,delta=c(1,1)),x$parameters)
  
  x <- weighted_j2r(trt.weight=1,delta=c(1.3,1.4))
  expect_equal(list(trt.weight=1,delta=c(1.3,1.4)),x$parameters)
  
})