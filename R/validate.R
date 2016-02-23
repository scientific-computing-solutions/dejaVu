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


#validates the covar argument to ensure it is not empty and
#does not include the treatment arm 
validatecovar <- function(covar){
  
  if(class(covar)!="formula") stop("Invalid covar is not of type formula")
  
  if(length(.getResponse(covar))!=0){
    stop("covar cannot have any variables on the left hand side.")
  }
  
  covariates <- attr(terms(covar),"term.labels")
  if(length(covariates)==0){
    stop("Empty covar argument!")
  }
  
  ans <- unlist(lapply(covariates,function(x){
    "arm" %in% unlist(strsplit(x,split=c(":")))
  }))
    
  if(any(ans)){
    stop("The covar argument cannot include interactions between treatment arm and covariates.")
  }

}

#from http://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
.getResponse <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response] 
}