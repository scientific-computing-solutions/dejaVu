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


#helper function for GetGamma_mu
.getPredict <- function(models,data){
  data$censored.time <- rep(1,nrow(data))
  data$arm <- factor(rep(0,nrow(data)))
  mu <- exp(predict(models[[1]],newdata=data,type="link"))
  data$arm <- factor(rep(1,nrow(data)))  
  c(mu,exp(predict(models[[2]],newdata=data,type="link")))  
}

#given model (or list of models) of from creating a SimFit object
#output a list of impute.parameters which contains
#the value of equal.dispersion
#a vector of gamma (1/dispersion) values, control then active arm
#a matrix of mu (mean of negative binomial) values, one row per subject
#first column control arm, second column active arm
GetGamma_mu <- function(model,data,equal.dispersion){
  if(equal.dispersion){
    mod.summary <- summary(model)
    gamma <- rep(mod.summary$theta,2)
    mu<- .getPredict(list(model,model),data)
  }
  else{
    gamma <- vapply(model,function(mod){summary(mod)$theta},FUN.VALUE = numeric(1)) 
    mu <-  .getPredict(model,data) 
  }
  #when this function is changed then don't forget to update SingleSimFit.object in Roxygen
  #and the user guide
  gamma_mu_function <- function(){list(gamma=gamma,mu=matrix(mu,ncol = 2,byrow=FALSE))}
  
  return(list(gamma_mu_function=gamma_mu_function,
              equal.dispersion=equal.dispersion))
}