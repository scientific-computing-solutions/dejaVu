#given a vector of subject arms (e.g. c(0,1,1,0)) and a vector
#for subject dispersions (e.g. 0.25 or c(0.28,0.25)) output a vector
#of subject specific dispersions where if dispersion is length 2 then
#subjects in arm 0 are given dispersions[1] and subjects in arm 1 are given
#dispersions[2]
getDispersions <- function(arm,dispersions){
  if(length(dispersions)==1){
    dispersions <- rep(dispersions,2)
  }
  ifelse(arm==0,dispersions[1],dispersions[2])
}

#Output a single subject's Poisson process rate given study.time,
#event.rate and dispersion parameter
GetSimRates <- function(study.time,event.rate,dispersion){
  if(dispersion==0){
    return(event.rate)  
  }
  lambda <- rgamma(n=1,shape=1/dispersion,scale=event.rate*study.time*dispersion)
  return(lambda/study.time)
}



# Given a single Poisson process rate and total study time
# returns a vector of event.times for a single subject
# if no events occur in the time [0,study.time] then numeric(0) 
# is returned
# This function does not take into account the requirement for the adjacent 
# exacerbations to be > 7 days apart
GetEventTimes <- function(rate,study.time){
  if(is.infinite(rate)){
    stop("An infinite rate cannot be used!")
  }
  event.times <- numeric(0)
  current.time <- rexp(1) / rate
  
  while(current.time <= study.time){
    event.times <- c(event.times,current.time)
    current.time <- current.time + rexp(1)/rate
  }
  
  return(event.times)
}


#If no formula is given to the fit function then 
#a default formula is used (different depending on whether the arms
#have the same dispersion parameters)
GetModelFormula <- function(equal.dispersion,covar){
  if(!is.logical(equal.dispersion) || length(equal.dispersion)>1){
    stop("Invalid equal.dispersion argument")
  }
  
  rhs <- if(equal.dispersion) "arm + offset(log(censored.time))" else "offset(log(censored.time))"
  
  if(!is.null(covar)) {
    validatecovar(covar)
    return(update.formula(covar,paste("observed.events~",rhs,"+.")))
  }  
  
  formula(paste("observed.events ~",rhs))
  
} 


#simple wrapper to create singleSim objects
.singleSimConstructor <- function(data, event.times, status, subject.rates,
                      dropout.mechanism, impute.mechanism,
                      study.time, event.rates, dispersions){

  retVal <- list(data=data, event.times=event.times, status=status, 
                 subject.rates=subject.rates,
                 dropout.mechanism=dropout.mechanism, 
                 impute.mechanism=impute.mechanism,
                 study.time=study.time, event.rates=event.rates,
                 dispersions=dispersions)
  class(retVal) <- "SingleSim"
  return(retVal)
}

#strip the dejaData object's data frame to leave only covariates
#(not subject Id, arm or rate)
remove.dejacols <- function(dejaData){
  dejaData$data[,dejaData$Id] <- NULL
  dejaData$data[,dejaData$arm] <- NULL
  if(!is.null(dejaData$rate)) dejaData$data[,dejaData$rate] <- NULL
  dejaData$data
}
