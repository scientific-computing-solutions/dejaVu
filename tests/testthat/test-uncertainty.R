context("uncertainty")

#Note tests with use.uncertainty=FALSE are included elsewhere

getSimFit <- function(equal.dispersion){
  sim <- SimulateComplete(study.time=365, 
                          number.subjects=100, 
                          event.rates=c(0.01,0.005),
                          dispersions=0.25) 
  
  Simfit(SimulateDropout(sim,drop.mechanism=ConstantRateDrop(rate=0.0025,var=1)),
         equal.dispersion=equal.dispersion)
}


test_that("apply.func",{
  
  ans <- apply.func(list(1,2,3),function(x,y)x*y,num.reps=3,y=2)
  expect_equal(ans,c(2,4,6,2,4,6,2,4,6))
  
  ans <- apply.func(list(2,5,3),function(x,y)x*y,num.reps=4,y=6)
  expect_equal(ans,c(12,30,18,12,30,18,12,30,18,12,30,18))
  
})


test_that("equal_disp_same_gamma",{
  set.seed(33143)
  sf <- getSimFit(TRUE)  
  
  x1 <- sf$genCoeff.function(TRUE)
  expect_equal(x1$gamma[1],x1$gamma[2])
  
})

test_that("not_equal_disp_diff_gamma",{
  set.seed(33143)
  sf <- getSimFit(FALSE)  
  
  x1 <- sf$genCoeff.function(TRUE)
  expect_false(x1$gamma[1]==x1$gamma[2])
  
})

test_that("same_mu_per_subjects",{
  set.seed(33143)
  sf <- getSimFit(TRUE)  
  
  x1 <- sf$genCoeff.function(TRUE)
  
  expect_true(all(x1$mu[,1]==x1$mu[1,1]))
  expect_true(all(x1$mu[,2]==x1$mu[1,2]))
})


test_that("is_stochastic",{
  
  set.seed(33143)
  sf <- getSimFit(TRUE)  
 
  x1 <- sf$genCoeff.function(TRUE)
  x2 <- sf$genCoeff.function(TRUE)
  
  expect_true(all(x1$gamma != x2$gamma))
  expect_false(x1$mu[1,1]==x2$mu[1,1])
})

test_that("takes_into_account_uncertainty",{
  
  set.seed(33143)
  sf <- getSimFit(TRUE)  
  
  x1 <- sf$genCoeff.function(TRUE)
  expect_warning(x2 <- sf$genCoeff.function(FALSE))
  
  expect_true(all(x1$gamma != x2$gamma))
  expect_false(x1$mu[1,1]==x2$mu[1,1])
})