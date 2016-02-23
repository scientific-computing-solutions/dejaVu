context("DropoutMechanism")

#Note no validation of GetDropTime function

test_that("CreateNewDropoutMechanism_invalid_args",{
  
  expect_error(CreateNewDropoutMechanism(type="MMAR",text="hello",GetDropTime=function(){}))
  expect_error(CreateNewDropoutMechanism(type=c("MAR","DE"),text="hello",GetDropTime=function(){}))
  expect_error(CreateNewDropoutMechanism(type="MAR",text=c("hello","x"),GetDropTime=function(){}))
  expect_error(CreateNewDropoutMechanism(type="MAR",text="hello",GetDropTime=45))
  
  expect_error(CreateNewDropoutMechanism(type="MAR",text="hello",GetDropTime=function(){},cols.needed=NULL))
  expect_error(CreateNewDropoutMechanism(type="MAR",text="hello",GetDropTime=function(){},parameters=c("hello")))
})

test_that("DummyDropout",{
  d <- CreateNewDropoutMechanism(type="MAR",text="hello",GetDropTime=function(){})
  expect_equal("DropoutMechanism",class(d))
  expect_equal("MAR",d$type)
  expect_equal("hello",d$text)
  expect_true(is.function(d$GetDropTime))
  
  d <- CreateNewDropoutMechanism(type="MNAR",text="hello",GetDropTime=function(){},cols.needed=c("a","b"),parameters=list(a=1,b=2))
  expect_equal("DropoutMechanism",class(d))
  expect_equal("MNAR",d$type)
  expect_equal("hello",d$text)
  expect_true(is.function(d$GetDropTime))
  expect_equal(c("a","b"),d$cols.needed)
  expect_equal(list(a=1,b=2),d$parameters)
})


test_that("ConstantRateDrop_validargs",{
  expect_error(ConstantRateDrop(rate=Inf))
  expect_error(ConstantRateDrop(rate=0))
  expect_error(ConstantRateDrop(rate=c(4,5)))
  expect_error(ConstantRateDrop(rate=5,var=-5))
  
  crd <- ConstantRateDrop(rate=0.1)
  expect_equal("MCAR",crd$type)
  expect_equal(list(rate=0.1,between.subject.var=0),crd$parameters)
  expect_equal("censored.time",crd$cols.needed)
  
  crd <- ConstantRateDrop(rate=0.1,var=4)
  expect_equal(list(rate=0.1,between.subject.var=4),crd$parameters)
})

test_that("LinearRateChangeDrop_validargs",{
  expect_error(LinearRateChangeDrop(starting.rate=-4,rate.change=0.5))
  expect_error(LinearRateChangeDrop(starting.rate=0,rate.change=0.5))
  expect_error(LinearRateChangeDrop(starting.rate=c(4,10),rate.change=1.5))
  expect_error(LinearRateChangeDrop(starting.rate=0.5,rate.change="dc"))
  expect_error(LinearRateChangeDrop(starting.rate=0.5,rate.change=0.005,var=-5))
  
  lrcd <- LinearRateChangeDrop(starting.rate=0.5,rate.change=-0.005)
  expect_equal("MAR",lrcd$type)
  expect_equal(list(starting.rate=0.5,rate.change.after.event=-0.005,between.subject.var=0),lrcd$parameters)
  expect_equal("censored.time",lrcd$cols.needed)
  
  lrcd <- LinearRateChangeDrop(starting.rate=0.5,rate.change=-0.005,var=4)
  expect_equal(list(starting.rate=0.5,rate.change.after.event=-0.005,between.subject.var=4),lrcd$parameters)
  
  
})


