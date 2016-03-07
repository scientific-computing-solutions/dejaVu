##' Create a weighted_j2r \code{ImputeMechanism} object
##' 
##' Missing counts for a subject in the active treatment arm will be imputed according to a point 
##' (determined by \code{trt.weight}) between the means of the placebo and treatment arms,
##' conditioned on the number of events. Missing counts for subjects in the 
##' placebo arm will be imputed according to the mean of the placebo arm, conditioned
##' on the subject's observed number of events. 
##' 
##' If \code{trt.weight} = 0 then imputation using this mechanism will follow 
##' the jump to reference (j2r) model whereby missing counts for subjects 
##' in both arms will be imputed according to the mean of the placebo arm 
##' conditioned on the subject's observed number of events
##' 
##' If \code{trt.weight} = 1 then imputation using this mechanism will follow 
##' the MAR model whereby missing counts for subjects 
##' in each arm will be imputed according to the event rate of subjects in its treatment group 
##' conditioned on the subject's observed number of events
##' 
##' See the User guide vignette for further details
##' 
##' @param trt.weight See details
##' @param delta If \code{trt.weight=1} then delta is a vector of length 2
##' (control.delta,treatment.delta) and the mean number of expected events for the imputed missing data is
##' multipled by the appropriate delta 
##' @return An \code{ImputeMechanism} object
##' @seealso \code{\link{ImputeMechanism.object}}
##' @export
weighted_j2r <- function(trt.weight,delta=c(1,1)){
  #first validate the arguments
  
  if(!.internal.is.finite.number(trt.weight) || trt.weight < 0 || trt.weight > 1){
    stop("Invalid argument trt.weight should be in [0,1]")
  }

  if(!is.numeric(delta) || any(is.na(delta)) || any(is.infinite(delta) || length(delta)!=2 || any(delta<=0))){
    stop("Invalid argument delta should be a vector of length 2 of positive numbers")
  }
  
  
  if(trt.weight!=1 && !all(delta==1)){
    stop("Invalid argument delta must be c(1,1) unless trt.weight=1")  
  }
  
  #define a helper function to be used in .internal.impute, using
  #the parameter trt.weight
  treatment.p.choice <- function(ps){
    return(ps[1]*(1-trt.weight)+ps[2]*(trt.weight))
  }
  
  #A function which takes SimFit object and outputs a list of 2 elements
  #1) newevent.times a list of vectors, the imputed event times for each subject (if subject has no new imputed
  #events then the vector should be numeric(0))
  #2) new.censored.times - the time at which subjects are censored in the imputed data set
  f <- function(fit){
    return(list(newevent.times=.internal.impute(fit,treatment.p.choice,delta),
                new.censored.times=rep(fit$singleSim$study.time,numberSubjects(fit))
    ))
    
  }
  
  CreateNewImputeMechanism(name="weighted_j2r",
                           cols.needed=c("censored.time","observed.events","arm"),
                           impute=f,
                           parameters=list(trt.weight=trt.weight,delta=delta))
}


.internal.impute <- function(fit,treatment.p.choice,delta){
  # performs the weighted_j2r method using the given SimFit object
  # a function which chooses the appropriate weighting on the treatment arm 
  # and the scaling factors delta
  # returns the imputed event times for each subjects as a list of vectors
  
  #@param fit is a SimFit object
  #@param treatment.p.choice is a function which given a vector of p's for a subject 
  #which returns the appropriately weighted value of p to be used in the imputation
  #@param delta argument in weighted_j2r function
  #@return returns a list of vectors, the imputed event times for each subject (if subject has no new imputed
  #events then the vector should be numeric(0))
  
  gamma_mu <- fit$impute.parameters$gamma_mu_function()
  
  #the data frame from the SimFit object  
  df <- fit$singleSim$data
  
  #take each subject in turn
  lapply(seq_along(nrow(df)), function(i){
    
    study.time <- fit$singleSim$study.time 
    time.left <- study.time - df$censored.time[i]
    if(time.left==0){return(numeric(0))} #subject was not censored
    
    if(df$arm[i]==0){
      p <- (gamma_mu$mu[i,1] * time.left)/(gamma_mu$gamma[1]+gamma_mu$mu[i,1]*study.time)
      gamma <- gamma_mu$gamma[1] + df$observed.events[i]
      delta.factor <- delta[1]
    }
    else{
      p <- (gamma_mu$mu[i,1]*time.left)/
           (gamma_mu$gamma[2] + gamma_mu$mu[i,2]*df$censored.time[i] + gamma_mu$mu[i,1]*time.left)
      p <- c(p, (gamma_mu$mu[i,2]*time.left)/(gamma_mu$gamma[2]+gamma_mu$mu[i,2]*study.time))
      p <- treatment.p.choice(p)
      gamma <- gamma_mu$gamma[2] + df$observed.events[i]
      delta.factor <- delta[2]
    }  
    
    u <- (p/(1-p))*gamma*delta.factor
    
    rate <- GetSimRates(time.left,event.rate=u/time.left,dispersion=1/gamma) 
    return(GetEventTimes(rate,time.left)+df$censored.time[i])  
    #note numeric(0)+x = numeric(0)
  })
  
}