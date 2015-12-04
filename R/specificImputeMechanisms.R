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
##' @param trt.weight See details
##' @param delta If \code{trt.weight=1} then delta is a vector of length 2
##' (control.delta,treatment.delta) and the mean number of expected events for the imputed missing data is
##' multipled by the appropriate delta 
##' @return An \code{ImputeMechanism} object
##' @seealso \code{\link{ImputeMechanism.object}}
##' @export
weighted_j2r <- function(trt.weight,delta=c(1,1)){
  
  if(!.internal.is.finite.number(trt.weight) || trt.weight < 0 || trt.weight > 1){
    stop("Invalid argument trt.weight should be in [0,1]")
  }

  if(!is.numeric(delta) || any(is.na(delta)) || any(is.infinite(delta) || length(delta)!=2 || any(delta<=0))){
    stop("Invalid argument delta should be a vector of length 2 of positive numbers")
  }
  
  
  if(trt.weight!=1 && !all(delta==1)){
    stop("Invalid argument delta must be c(1,1) unless trt.weight=1")  
  }
  
  treatment.p.choice <- function(ps){
    return(ps[1]*(1-trt.weight)+ps[2]*(trt.weight))
  }
  
  
  f <- function(fit){
    return(list(newevent.times=.internal.impute(fit,treatment.p.choice,delta),
                new.censored.times=rep(fit$singleSim$study.time,numberSubjects(fit))
    ))
    
  }
  
  retVal <- list(name="weighted_j2r",
                 cols.needed=c("censored.time","observed.events","arm"),
                 impute=f,
                 parameters=list(trt.weight=trt.weight,delta=delta))
  class(retVal) <- "ImputeMechanism"
  return(retVal)
}

# performs the weighted_j2r method using the given SimFit object
# a function which chooses the appropriate weighting on the treatment arm 
# and the scaling factors delta
# returns the imputed event times for each subjects as a list of vectors
.internal.impute <- function(fit,treatment.p.choice,delta){
  
  df <- fit$singleSim$data
  
  lapply(1:nrow(df), function(i){
    study.time <- fit$singleSim$study.time 
    time.left <- study.time - df$censored.time[i]
    if(time.left==0){return(numeric(0))}
    
    if(df$arm[i]==0){
      p <- (fit$impute.parameters$p[1] * time.left)/(fit$impute.parameters$gamma[1]+fit$impute.parameters$p[1]*study.time)
      gamma <- fit$impute.parameters$gamma[1] + df$observed.events[i]
      delta.factor <- delta[1]
    }
    else{
      p <- (fit$impute.parameters$p[1]*time.left)/
           (fit$impute.parameters$gamma[2] + fit$impute.parameters$p[2]*df$censored.time[i] + fit$impute.parameters$p[1]*time.left)
      p <- c(p, (fit$impute.parameters$p[2]*time.left)/(fit$impute.parameters$gamma[2]+fit$impute.parameters$p[2]*study.time))
      p <- treatment.p.choice(p)
      gamma <- fit$impute.parameters$gamma[2] + df$observed.events[i]
      delta.factor <- delta[2]
    }  
    
    u <- (p/(1-p))*gamma*delta.factor
    
    rate <- GetSimRates(time.left,number.subject=1,event.rate=u/time.left,dispersion=1/gamma) 
    return(GetEventTimes(rate,time.left)+df$censored.time[i])  
    #note numeric(0)+x = numeric(0)
  })
  
}