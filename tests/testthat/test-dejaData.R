context("DejaData")

test_that("defaultDejaData",{
  expect_error(defaultDejaData(50,"sdf"))
  expect_error(defaultDejaData(50,0))
  expect_error(defaultDejaData(50,c(4,-3)))
  expect_error(defaultDejaData(50,c(0.5,0.6,0.7)))
  expect_error(defaultDejaData(0,c(0.5,0.6)))
  expect_error(defaultDejaData(c(10,40,50),c(0.5,0.6)))
  
  dd <- defaultDejaData(25,c(0.75,0.34))
  expect_equal("arm",dd$arm)
  expect_equal("Id",dd$Id)
  expect_equal("rate",dd$rate)
  
  expect_equal(c(rep(0,25),rep(1,25)),dd$data$arm)
  expect_equal(c(rep(0.75,25),rep(0.34,25)),dd$data$rate)
  expect_equal(1:50,dd$data$Id)

  dd <- defaultDejaData(c(10,15),c(0.75,0.34))
  expect_equal(c(rep(0,10),rep(1,15)),dd$data$arm)
  expect_equal(c(rep(0.75,10),rep(0.34,15)),dd$data$rate)
  
  dd <- defaultDejaData(c(10,15),0.75)
  expect_equal(c(rep(0,10),rep(1,15)),dd$data$arm)
  expect_equal(rep(0.75,25),dd$data$rate)

})

test_that("MakeDejaData_valid",{
  df <- data.frame(Id=1:4,myarm=c(0,1,1,0),myrate=c(1,2,3,4),X=c(0,1,1,1))
  dd <- MakeDejaData(data=df,arm="myarm",rate="myrate",Id="Id")

  expect_equal("DejaData",class(dd))
  expect_equal(df,dd$data)
  expect_equal("myarm",dd$arm)
  expect_equal("myrate",dd$rate)
  expect_equal("Id",dd$Id)
  
  #OK with rate =NULL
  dd <- MakeDejaData(data=df,arm="myarm",Id="Id")
  expect_equal("DejaData",class(dd))
  expect_equal(df,dd$data)
  expect_true(is.null(dd$rate))
})

test_that("MakeDejaData_invalid",{
  
  df <- data.frame(Id=1:4,myarm=c(0,1,1,0),myrate=c(1,2,3,4),X=c(0,1,1,1))
  
  #not a data frame
  expect_error(MakeDejaData(data="hello",arm="myarm",rate="myrate",Id="Id"))
  
  #invalid column names
  expect_error(MakeDejaData(data=df,arm="arm",rate="myrate",Id="Id"))
  expect_error(MakeDejaData(data=df,arm="myarm",rate="m",Id="Id"))
  expect_error(MakeDejaData(data=df,arm="myarm",rate="myrate",Id="ID"))
  
  #Ids not unique
  df$Id[1] <- 3
  expect_error(MakeDejaData(data=df,arm="myarm",rate="myrate",Id="Id"))
  df$Id[1] <- 1
  
  #arm doesn't include both 0 and 1
  df$myarm <- rep(0,4)
  expect_error(MakeDejaData(data=df,arm="myarm",rate="myrate",Id="Id"))
  df$myarm <- rep(1,4)
  expect_error(MakeDejaData(data=df,arm="myarm",rate="myrate",Id="Id"))
  df$myarm <- c(0,1,1,0)
  
  #arm includes something other than 0 and 1
  df$myarm[1] <- 4
  expect_error(MakeDejaData(data=df,arm="myarm",rate="myrate",Id="Id"))
  df$myarm[1] <- 0
  
  #rate is not positive
  df$myrate[1] <- -5
  expect_error(MakeDejaData(data=df,arm="myarm",rate="myrate",Id="Id"))
  
  
})

