context("CopyReference")

test_that("copyReference",{
  skip_on_cran()
  expect_equal({
    set.seed(601231)
    n <- 10000
    z <- c(rep(0,n/2),rep(1,n/2))
    
    k <- 4
    gamma <- rgamma(n,shape=k,scale=1/k)
    rate <- gamma*exp(1*z)
    
    #simulate random follow-up time for subject if they remain in study
    fullFupTime <- 1+0.1*runif(n)
    
    #simulate time at which they might drop out, between 0 and 1
    dropOutTime <- runif(n)
    
    y1 <- rpois(n, lambda=rate*dropOutTime)
    #generate count for second period assuming they remain on randomised treatment
    y2 <- rpois(n, lambda=rate*(fullFupTime-dropOutTime))
    #generate another outcome that occurs if patient is on placebo in second period,
    #using the copy reference assumption
    upsi1 <- dropOutTime*exp(0)
    upsi2 <- (fullFupTime-dropOutTime)*exp(0)
    y2placebo <- rnbinom(n, size=(k+y1), prob=((k+upsi1)/(k+upsi1+upsi2)))
    
    #model dropout process, conditional on y1
    do_xb <- -1+1*log(y1+1)
    do_pr <- exp(do_xb)/(1+exp(do_xb))
    do <- 1*(runif(n)<do_pr)
    
    y2[(do==1) & (z==1)] <- y2placebo[(do==1) & (z==1)]
    
    y <- y1+y2*(1-do)
    fupTime <- fullFupTime
    fupTime[do==1] <- dropOutTime[do==1]
    
    #analyse full data - the treatment effect from this is what the copy reference
    #MI analysis is targeting
    fullDataEst <- MASS::glm.nb(y1+y2~z+offset(log(fullFupTime)))$coef[2]
    
    #save data to file for comparison with estimate from SAS
    #write.csv(data.frame(z=z,y=y,time=fupTime), file="copyRefTest.csv",row.names=FALSE)
    #running using James Roger's macros (from www.missingdata.org.uk), we obtain
    #a very close point estimate and 95% CI to what is obtained below using dejavu
    
    #perform MI under the copy reference assumption
    covar.df <- data.frame(id=1:n, z=z)
    dejaData <- MakeDejaData(covar.df, arm="z", Id="id")
    
    #load data using ImportSim
    obsData <- ImportSim(dejaData, event.times=expandEventCount(count=y, time=fupTime), status="dropout",
                         study.time=1, censored.time=fupTime, allow.beyond.study=T)
    obsFit <- Simfit(obsData,equal.dispersion=TRUE)
    
    imputed.data.sets <- Impute(fit = obsFit,impute.mechanism = copy_reference(),
                                N=10)
    #analyse imputed datasets
    fitted <- Simfit(imputed.data.sets,family="negbin")
    crMIEst <- log(summary(fitted)$treatment.effect)
  
    #compare full data estimates with MI estimate
    as.logical(abs(fullDataEst-crMIEst)<0.01)
  }, TRUE)
})