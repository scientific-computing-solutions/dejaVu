context("DropoutMechanism")

#Note no validation onGetDropTime function

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