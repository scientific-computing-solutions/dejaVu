# Validates the user arguments to SimulateComplete function
# for arguments see SimulateComplete function
# if invalid argument an exception is thrown
ValidateSimCompleteArgs <- function(study.time,number.subjects,event.rates,dispersions){
  
  if(!.internal.is.finite.number(study.time) || study.time < 0){
    stop("Invalid study.time argument it must be a single positive finite numeric value")
  }
  
  lapply(c(number.subjects,event.rates,dispersions),function(x){
    if(!is.numeric(x) || x<0|| is.na(x) || is.infinite(x)){
      stop("Invalid argument: number.subjects, event.rates and dispersions must be non-negative vectors of length at most 2")
    }
  })
  
  if(length(number.subjects) > 2 ||length(event.rates) > 2 || length(dispersions) > 2 ){
    stop("Invalid argument: number.subjects, event.rates and dispersions must be non-negative vectors of length at most 2")
  }
  
  if(!all(.internal.is.wholenumber(number.subjects))){
    stop("Invalid argument: there must be a integer number of subjects")
  }
  
  if(any(event.rates==0)|| any(number.subjects==0)){
    stop("Invalid argument: number.subjects and event.rates cannot be zero")
  }
  
  
}

# Return a vector of length number.subject of rates for the Poission process
# 1 for each subject of a given arm. The arguments here are single values not
# numeric vectors
GetSimRates <- function(study.time,number.subject,event.rate,dispersion){
  if(dispersion==0){
    return(rep(event.rate,number.subject))  
  }
  lambda <- rgamma(n=number.subject,shape=1/dispersion,scale=event.rate*study.time*dispersion)
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


# validate the SimulateDropout function
# for argumenets see SimulateDropout function
# if invalid argument an exception is thrown
validateSimulateDropout <- function(simComplete,drop.mechanism){
  if(class(simComplete)!="SingleSim"){
    stop("Invalid argument: simComplete must be of class SingleSim")
  }
  
  if(class(drop.mechanism)!="DropoutMechanism"){
    stop("Invalid argument: drop.mechanism must be of class DropoutMechanism")
  }
  
  if(simComplete$status!="complete"){
    stop("Invalid SingleSim argument: simComplete$status != 'complete'")
  }
  
  if(!all(drop.mechanism$cols.needed %in% colnames(simComplete$data))){
    stop("This dropout mechanism requires ",paste(drop.mechanism$cols.needed,collapse=", "),
         "as column names in the simulated data frame")
  }
  
}

#If no formula is given to the fit function then 
#a default formula is used (different depending on whether the arms
#have the same dispersion parameters)
GetDefaultFormula <- function(equal.dispersion){
  if(!is.logical(equal.dispersion) || length(equal.dispersion)>1){
    stop("Invalid equal.dispersion argument")
  }
  
  if(equal.dispersion){
    return(formula("observed.events ~ arm + offset(log(censored.time))"))
  }
  
  formula("observed.events ~ offset(log(censored.time))")
  
} 

#validate the arguments for the SimFit.SingleSim function
ValidateSimFitArguments <- function(family,equal.dispersion){
  
  if(!is.logical(equal.dispersion) || length(equal.dispersion)>1){
    stop("Invalid argument: equal.dispersion")
  }

  allowed.models <- c("negbin","poisson","quasipoisson")
  if(!family %in% allowed.models){
    stop("Invalid argument: family must be one of",paste(allowed.models,sep=", "))
  }
  
  if(family!="negbin" && !equal.dispersion){
    stop("equal.dispersion argument cannot be used unless negative binomial model is to be fit")
  }
  
}
