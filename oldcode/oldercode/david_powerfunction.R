##########################
## MODIFICATION DAVID RUAU
##########################

## libraries
library(MASS)
library(plyr)

## libraries to paralellise the loops
library(foreach)
library(doMC)

# detecting the number of core on the current computer  
ncore = detectCores()
registerDoMC(cores = ncore)

## using object to store your results
## This is an example I did not had time to implement the 
## use of this code int he for loop

##' Represent the result of running a simulation
##' @se Standard error
##' @se2 Standard error two
##' @rratio REPLACE ME
##' @pval pvalue
##' @theta REPLACE ME
##' @warn warning

setClass("distribution",
       representation(se="numeric",
                      se2="numeric",
                      rratio="numeric",
                      pval="numeric",
                      theta="numeric",
                      warn="numeric"))
  
complete <- function(n, rate.placebo, rr, T, dispersion.placebo, 
dispersion.treatment, N.rep, alpha){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  
  ## this is how you would create your objects.
  ## You create a object of class 'distribution' containing all
  ## the fields necessary
  poisson <- new("distribution")
  quasipoisson <- new("distribution")
  nbd <- new("distribution")
  
  err <- rep(NA,N.rep)
		
  nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
  nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)

  dropout_time.placebo <- rep(T,n)
  dropout_time.treatment <- rep(T,n)
    
  nothing <- foreach(i=1:N.rep, .combine=c) %dopar% {
		
	## this number does not change. it is a constant number.
  	dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
  	dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
  
  	nb_incomplete.placebo   <- nb_incomplete(n, nb_complete.placebo,   dropout_time.placebo)
  	nb_incomplete.treatment <- nb_incomplete(n, nb_complete.treatment, dropout_time.treatment)

    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    poisson.rratio[i] <- exp(poisson.summary$coeffi[2,1])
    poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
    poisson.se2[i] <- poisson.summary$coeffi[2,2]
    poisson.pval[i] <- poisson.summary$coef[2,4] #p-value
    
    quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    quasipoisson.rratio[i] <- exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]
    quasipoisson.pval[i] <- quasipoisson.summary$coef[2,4] 
    
    nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    nbd.warn[i] <- (!is.null(nbd.summary$th.warn))
    nbd.rratio[i] <- exp(nbd.summary$coeffi[2,1])
    nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])
    nbd.se2[i] <- nbd.summary$coeffi[2,2]
    nbd.pval[i] <- nbd.summary$coef[2,4] 
    nbd.theta[i] <- nbd.summary$theta
}

