# Validates the user arguments to create the DejaData object
# within the SimulateComplete function
# for arguments see SimulateComplete function
# if invalid argument an exception is thrown
validateDefaultDejaDataArgs <- function(number.subjects,event.rates){
  
  lapply(c(number.subjects,event.rates),function(x){
    if(is.null(x) || !is.numeric(x) || x<0|| is.na(x) || is.infinite(x)){
      stop("Invalid argument: number.subjects, event.rates and dispersions must be non-negative vectors of length at most 2")
    }
  })
  
  if(length(number.subjects) > 2 ||length(event.rates) > 2 ){
    stop("Invalid argument: number.subjects, event.rates must be non-negative vectors of length at most 2")
  }
  
  if(!all(.internal.is.wholenumber(number.subjects))){
    stop("Invalid argument: there must be a integer number of subjects")
  }
  
  if(any(event.rates==0)|| any(number.subjects==0)){
    stop("Invalid argument: number.subjects and event.rates cannot be zero")
  }
}

#validate the arguments to simcomplete which are not needed to create
#a default dejaData argument. If invalid argument an exception is thrown
ValidateSimCompleteArgs <- function(dejaData,study.time,dispersions){
  if(!.internal.is.finite.number(study.time) || study.time < 0){
    stop("Invalid study.time argument it must be a single positive finite numeric value")
  }
  
  if(class(dejaData) != "DejaData"){
    stop("Invalid dejaData argument")
  }
  
  if(is.null(dejaData$rate)){
    stop("Cannot use a DejaData object without a rate column when simulating")
  }
  
  if(length(dispersions) > 2){
    stop("Invalid argument: dispersions must be non-negative vectors of length at most 2")
  }
  
  lapply(dispersions,function(x){
    if(!is.numeric(x) || x<0|| is.na(x) || is.infinite(x)){
      stop("Invalid argument: dispersions must be non-negative vectors of length at most 2")
    }
  })
  
}

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

  if(!is.character(family)){
    stop("Invalid argument: family must be a character")
  }
  
  allowed.models <- c("negbin","poisson","quasipoisson")
  if(!family %in% allowed.models){
    stop("Invalid argument: family must be one of",paste(allowed.models,sep=", "))
  }
  
  if(family!="negbin" && !equal.dispersion){
    stop("equal.dispersion argument cannot be used unless negative binomial model is to be fit")
  }
  
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