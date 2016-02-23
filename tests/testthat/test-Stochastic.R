context("stochastic_tests")

test_that("negative_binomial",{
  require(stats)

  set.seed(104)
  CreateSamples <- function(p,gamma,N,Ti=100){
    dispersion <- 1/gamma
    event.rate <- p*gamma/((1-p)*Ti)
    SimulateComplete(study.time = Ti,number.subjects=N, event.rates=event.rate,dispersions=dispersion)$data$observed.events
  }
  
  #only stochastically true so be careful
  expect_warning(expect_true(ks.test(CreateSamples(p=0.5,gamma=1,N=10000),"pnbinom",size=1,prob=0.5)$p.value<0.05))
  expect_warning(expect_true(ks.test(CreateSamples(p=0.25,gamma=1.3,N=10000,Ti=50),"pnbinom",size=1.3,prob=0.25)$p.value<0.05))
  expect_warning(expect_true(ks.test(CreateSamples(p=0.15,gamma=0.8,N=10000),"pnbinom",size=0.8,prob=0.15)$p.value<0.05))
  
})

test_that("poisson",{
  require(stats)
  
  set.seed(1304)
  CreateSamples <- function(rate,N,Ti=100){
    SimulateComplete(study.time = Ti,number.subjects=N, event.rates=rate,dispersions=0)$data$observed.events
  }
  
  expect_warning(expect_true(ks.test(CreateSamples(rate=0.025,N=10000,Ti=100),"ppois",lambda=100*0.025)$p.value<0.05))
  
})