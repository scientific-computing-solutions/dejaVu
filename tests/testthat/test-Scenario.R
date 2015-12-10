context("scenario")

test_that("CreateScenario",{
  
  #try with dummy classes 
  
  x <- list(list(2,3,4))
  expect_error(CreateScenario(x))
  
  class(x[[1]]) <- "summary.SingleSimFit"
  expect_error(CreateScenario(x,description=c("s","£sd")))
  
  x[[2]] <- list(3,4,5)
  class(x[[2]]) <- "summary.ImputeSimFit" 
  expect_error(CreateScenario(x))
  
  class(x[[1]]) <- "summary.ImputeSimFit" 
  
  cs <- CreateScenario(x,description="meh")
  expect_equal("meh",cs$description)
  expect_equal(x,cs$summaries)
  
  class(x[[1]]) <- "summary.SingleSimFit"
  class(x[[2]]) <- "summary.SingleSimFit"
  cs <- CreateScenario(x,description="meh")
  expect_equal("meh",cs$description)
  expect_equal(x,cs$summaries)
})


example.scenario <- function(){ 
  
  #simulate a complete data set
  sim <- SimulateComplete(study.time=365,number.subjects=100,
                          event.rates=c(0.01,0.005),dispersions=0.25)
  
  #take the simulated data set and apply an MCAR dropout mechanism...                         
  sim.with.MCAR.dropout <- SimulateDropout(sim,
                                           drop.mechanism=ConstantRateDrop(rate=0.0025)) 
  
  #fit a Negative Binomial model 
  with.MCAR.fit <- Simfit(sim.with.MCAR.dropout,equal.dispersion=TRUE)
  
  #we can impute a set of 10 sets following the j2r mechanism using the fit
  impute.data.sets <- Impute(with.MCAR.fit,impute.mechanism = weighted_j2r(trt.weight=0),N=10)
  
  #we can then fit models to the entire imputed data set
  fit.imputed.set <- Simfit(impute.data.sets)
  
  #output the summary values
  return(list(MI=summary(fit.imputed.set), #for MI
              complete=summary(Simfit(sim,family="poisson")))) #for complete data set
}

test_that("create_scenario",{
  
  set.seed(10)
  
  ans <- replicate(2,example.scenario(),simplify = FALSE)
  s1 <- CreateScenario(list(ans[[1]]$complete,ans[[2]]$complete))
  s2 <- CreateScenario(list(ans[[2]]$MI,ans[[2]]$MI),description="hello")
  
  expect_equal("Scenario",class(s1))
  expect_equal("Scenario",class(s2))
  
  expect_equal("",s1$description)
  expect_equal("hello",s2$description)
  
  expect_equal(ans[[1]]$complete,s1$summaries[[1]])
  expect_equal(ans[[2]]$MI,s2$summaries[[2]])
  
})


test_that("as.data.frame.scenario",{
  #integration test
  set.seed(48)
  
  ans <- replicate(2,example.scenario(),simplify = FALSE)
  s1 <- CreateScenario(list(ans[[1]]$complete,ans[[2]]$complete))
  s2 <- CreateScenario(list(ans[[2]]$MI,ans[[2]]$MI),description="hello")
  
  expect_error(as.data.frame(s1,use.adjusted.pval = TRUE))
  
  df1 <- as.data.frame(s1)
  df2 <- as.data.frame(s2)
  df3 <- as.data.frame(s2,use.adjusted.pval = TRUE)
  
  expect_equal(2,nrow(df1))
  expect_equal(1:2,df1$replica)
  expect_equal(s1$summaries[[1]]$dropout[1],df1$dropout.control[1] )
  expect_equal(c(0,0),df1$dropout.active )
  
  expect_equal(c(198,198),df1$df)
  expect_equal(s1$summaries[[1]]$treatment.effect,df1$treatment.effect[1])
  
  expect_equal(df2$df[1],s2$summaries[[1]]$df)
  expect_equal(df3$df[1],s2$summaries[[1]]$adjusted.df)
  expect_equal(df2$pval[2],s2$summaries[[2]]$pval)
  expect_equal(df3$pval[2],s2$summaries[[2]]$adjusted.pval)
  expect_equal(df2$dispersion[1],s2$summaries[[1]]$dispersion)
  
  expect_equal(df2$dropout.active[1],s2$summaries[[1]]$dropout[2])

})

test_that("summary.Scenario.args",{
  set.seed(48)
  
  ans <- replicate(2,example.scenario(),simplify = FALSE)
  s1 <- CreateScenario(list(ans[[1]]$complete,ans[[2]]$complete))
  s2 <- CreateScenario(list(ans[[2]]$MI,ans[[2]]$MI),description="hello")
  
  expect_error(summary(s1,use.adjusted.pval = TRUE))
  expect_error(summary(s2,alpha=-5))
  expect_error(summary(s2,alpha=0))
  expect_error(summary(s2,alpha=c(0.1,0.3)))
  
})

test_that("internal.summary.Scenario",{
  
  my.df <- data.frame(treatment.effect=c(5,6,7),
                      se=c(0.6,0.9,1.4),
                      pval=c(0.6,0.06,0.01),
                      dropout.control=c(1,2,3),
                      dropout.active=c(4,5,6))
  
  .internal.summary.Scenario(my.df,alpha=0.05)
  
})
