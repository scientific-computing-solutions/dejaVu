# Validates the user arguments to SimulateComplete function
# for arguments see SimulateComplete function
# if invalid argument an exception is thrown
ValidateSimCompleteArgs <- function(study.time,number.subjects,event.rates,dispersions){
  
  if(!is.numeric(study.time) || length(study.time) > 1 || study.time < 0 || is.infinite(study.time) || is.na(study.time)){
    stop("Invalid study.time argument it must be a single positive finite numeric value")
  }
  
  lapply(c(number.subjects,event.rates,dispersions),function(x){
    if(length(x)>2 || !is.numeric(x) || any(x<0)|| any(is.na(x)) || any(is.infinite(x))){
      stop("Invalid argument: number.subjects, event.rates and dispersions must be non-negative vectors of length at most 2")
    }
  })
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if(!all(is.wholenumber(number.subjects))){
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
# This functiondoes not take into account the requirement for the adjacent 
# exacerbations to be > 7 days apart
GetEventTimes <- function(rate,study.time){
  event.times <- numeric(0)
  current.time <- rexp(1,rate)
  
  while(current.time <= study.time){
    event.times <- c(event.times,current.time)
    current.time <- current.time + rexp(1,rate)
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
  
  ValidColumns(drop.mechanism,colnames(simComplete$data))
  
}
