##' Simulate a complete data set
##' 
##' Simulate a complete data set of a recurrent event
##' clinical trial without dropouts using a negative binomial model
##' with given rates and dispersion parameters
##' 
##' Each subject's events are described by a Poisson process with a subject specific rate given by
##' \code{lambda/study.time} where \code{study.time} is the study follow up period and \code{lambda}
##' has a gamma distribution with \code{shape=1/dispersion} and \code{scale=dispersion*event.rate*study.time} 
##' 
##' Different dispersions, event.rates and number of subjects can be specified for both arms of the trial 
##' 
##' @param study.time The study follow up period
##' @param dejaData If not NULL this should contain a \code{DejaData} object. If this is used then
##' number.subjects and event.rates arguments are ignored
##' @param number.subjects The number of subjects, if a vector \code{c(a,b)} then
##' \code{a} subjects on the control arm and \code{b} subjects on the active arm. If 
##' \code{number.subjects} is a single number then both arms have the given number of subjects.
##' @param event.rates The rate parameter(s) for the negative binomial model (if single parameter then it is used for both arms)
##' @param dispersions The dispersion parameter(s) for the negative binomial model (if single parameter then it is used for both arms)
##' @return A \code{SingleSim} object with \code{status='complete'} 
##' @seealso \code{\link{SingleSim.object}} 
##' @examples
##' sim <- SimulateComplete(study.time=365,number.subjects=50,
##'                         event.rates=c(0.01,0.005),dispersions=0.25)
##' sim2 <- SimulateComplete(study.time=365,number.subjects=c(50,75),
##'                         event.rates=c(0.01,0.005),dispersions=c(0,0.25))
##' @export
SimulateComplete <- function(study.time,dejaData=NULL,number.subjects=NULL,event.rates=NULL,dispersions){
  
  if(is.null(dejaData)){
    dejaData <- defaultDejaData(number.subjects,event.rates)
  }
  
  ValidateSimCompleteArgs(dejaData,study.time,dispersions)
  
  subject.rates <- unlist(mapply(GetSimRates,study.time=study.time,
                          event.rate=dejaData$data[,dejaData$rate],
                          dispersion=getDispersions(dejaData$data[,dejaData$arm],dispersions),
                          SIMPLIFY=FALSE))

  event.times <- lapply(subject.rates,GetEventTimes,study.time=study.time)
  
  events <- vapply(event.times,length,FUN.VALUE=numeric(1))
  
  data <- data.frame(Id=dejaData$data[,dejaData$Id],
                     arm=as.factor(dejaData$data[,dejaData$arm]),   
                     censored.time=rep(study.time,nrow(dejaData$data)),
                     observed.events=events,
                     actual.events=events)
  
  if(ncol(dejaData$data)>3){
    dejaData$data[,dejaData$Id] <- NULL
    dejaData$data[,dejaData$arm] <- NULL
    dejaData$data[,dejaData$rate] <- NULL
    data <- cbind(data,dejaData$data)
  }
  
  
  retVal <- list(data=data,
                 event.times=event.times,
                 status="complete",
                 subject.rates=subject.rates, #needed? 
                 dropout.mechanism=NULL,
                 impute.mechanism=NULL,
                 study.time=study.time,
                 event.rates=event.rates,
                 dispersions=dispersions)
                 
  class(retVal) <- "SingleSim"
  return(retVal)
}


##' @export
print.SingleSim <- function(x,...){
  
  cat(x$status,"dataset\n")
  cat("Study follow up period:",x$study.time,"\n")
  cat("Negative binomial event rates:",x$event.rates,"\n")
  cat("Negative binomial dispersion:",x$dispersions,"\n")
  if(x$status!="complete"){
    print(x$dropout.mechanism)
    cat("\n")
  }
  if(x$status=="imputed"){
    print(x$impute.mechanism)
    cat("\n")
  }
  cat("Data:\n")
  str(x$data)
}


##' SingleSim Object
##' 
##' A class containing the data for a single simulation. Depending on the value
##' of \code{status}, this may be a complete data set, a set including subject dropouts
##' or a data set after multiple imputation  
##'   
##' \code{print.SingleSim} and \code{summary.SingleSim} methods are defined.
##'
##' @section Structure: The above components must be included in
##' a SingleSim Object
##' @param data The data frame, one row per subjecy containing (at least) the following columns
##' Id, arm, censored.time, observed.events and actual.events
##' @param event.times A list of event times. event.times[[1]] is a list of event times for subject with Id 1
##' The length of event.times[[1]] = the number of observed events of subject with Id 1 
##' @param status Either "complete", "dropout" or "imputed" denoting the status of the data set.
##' @param subject.rates A vector of the specific rates used for the Poisson process for subjects when generating the data
##' @param dropout.mechanism If status is not "complete" then this contains the \code{DropoutMechanism} object
##' used to perform the subject dropout. See \code{\link{DropoutMechanism.object}}.
##' @param impute.mechanism If the status is "imputed" then this contains the \code{ImputeMechanism} object
##' used to perform the imputation. See \code{\link{ImputeMechanism.object}}
##' @param study.time The study follow up period (see \code{SimulateComplete})
##' @param event.rates The control/active event rates (see \code{SimulateComplete}), if data set was 
##' generated without using these (e.g. the dejaData argument was used) then this is set to NULL 
##' @param dispersions  The control/active dispersion rates (see \code{SimulateComplete})
##' @name SingleSim.object
##' @aliases print.SingleSim summary.SingleSim
NULL


##' Simulate subject dropout
##' 
##' This function takes a complete recurrent event data set
##' and drop out mechanism and creates a data set set with
##' dropout
##'  
##' @param simComplete A \code{SingleSim} object (with \code{status="complete"})
##' @param drop.mechanism A \code{DropoutMechanism} object
##' @return A \code{SingleSim} object with \code{status='dropout'}
##' @examples 
##' sim <- SimulateComplete(study.time=365,number.subjects=50,
##'                         event.rates=c(0.01,0.005),dispersions=0.25)
##'                         
##' sim.with.MCAR.dropout <- SimulateDropout(sim,
##'                      drop.mechanism=ConstantRateDrop(rate=0.0025)) 
##' sim.with.MAR.dropout <- SimulateDropout(sim,
##'                      drop.mechanism=LinearRateChangeDrop(starting.rate=0.0025,rate.change=0.0005))                       
##'                                                                    
##' @export 
SimulateDropout <- function(simComplete,drop.mechanism){
  validateSimulateDropout(simComplete,drop.mechanism)
  
  censored.time <- vapply(seq_along(simComplete$event.times),function(i){
    drop.mechanism$GetDropTime(event.times=simComplete$event.times[[i]],data=simComplete$data[i,])
  },FUN.VALUE=numeric(1))
  
  simComplete$event.times <- mapply(function(censor.time,event.times){ 
                          return(event.times[event.times<=censor.time]) },
                        censor.time=censored.time,event.times=simComplete$event.times,SIMPLIFY = FALSE)
  
  simComplete$data$observed.events <- vapply(simComplete$event.times,length,FUN.VALUE=numeric(1))
  simComplete$data$censored.time <- censored.time
  simComplete$status <- "dropout"
  simComplete$dropout.mechanism <- drop.mechanism
  return(simComplete)
}


##' @export
numberSubjects.SingleSim <- function(x){
  return(nrow(x$data))  
}

##' @export
Simfit.SingleSim <- function(x,family="negbin",equal.dispersion=TRUE,formula=GetDefaultFormula(equal.dispersion=equal.dispersion),...){
  
  ValidateSimFitArguments(family,equal.dispersion) #No formula validation yet
  
  data <- x$data[x$data$censored.time>0,]
  gamma <- NULL
  p <- NULL
 
  if(equal.dispersion && family=="negbin"){
    model <- glm.nb(formula=formula,data=data,...)
    mod.summary <- summary(model)
    gamma <- rep(mod.summary$theta,2)
    p <- exp(mod.summary$coeffi[1,1])*c(1,exp(mod.summary$coeffi[2,1]))
  } else if(equal.dispersion){
    model <- glm(formula=formula,data=data,family=do.call(family,args=list()),...)
  }
  else{ #only in negbinomial case
    model <- lapply(0:1,function(i){glm.nb(formula=formula,data=data[data$arm==i,],...)})
    gamma <- vapply(model,function(mod){summary(mod)$theta},FUN.VALUE = numeric(1)) 
    p <-     vapply(model,function(mod){exp(summary(mod)$coeffi[1,1])},FUN.VALUE = numeric(1))
  }
  
  if(family=="negbin"){
    impute.parameters <- list(gamma=gamma,p=p,equal.dispersion=equal.dispersion)
    class(impute.parameters) <- "ImputeParameters"
  }
  else{
    impute.parameters <- list(equal.dispersion=equal.dispersion)
  }
    
  retVal <- list(singleSim=x,
                 model=model,
                 impute.parameters=impute.parameters)
  class(retVal) <- "SingleSimFit"
  return(retVal)
}


##' @export
summary.SingleSim <- function(object,...){
  
  .extract <- function(f){
    vapply(unique(object$data$arm),
         function(arm)f(object$data[object$data$arm==arm,]),FUN.VALUE = numeric(1))
  }
  
  total.events <- .extract(function(y){sum(y$observed.events)})
  time.at.risk <- .extract(function(y){sum(y$censored.time)})
  
  retVal <- list(status=object$status,
                 study.time=object$study.time,
                 number.subjects=.extract(nrow),
                 number.dropouts=.extract(function(y){sum(y$censored.time!=object$study.time)}),
                 total.events=total.events,
                 time.at.risk=time.at.risk,
                 empirical.rates=total.events/time.at.risk)
  class(retVal) <- "summary.SingleSim"
  return(retVal)
}


##' @export
print.summary.SingleSim <- function(x,...){
  cat(x$status,"dataset\n")
  cat("Study follow up period:",x$study.time,"\n")
  cat("Subjects (per arm):",x$number.subjects,fill = TRUE)
  cat("Subject dropouts (per arm):",x$number.dropouts,fill = TRUE)
  cat("Number of events (per arm):",x$total.events,fill=TRUE)
  cat("Total Time at risk (per arm):",x$time.at.risk,fill=TRUE)
  cat("Empirical event rates (per arm):",x$empirical.rates,fill=TRUE)
}


##' summary.SingleSim object
##' 
##' The object returned when calling the summary function on a \code{SingleSim} object
##' 
##' @param status The status of the SingleSim object
##' @param study.time The study.time from the SingleSim object
##' @param number.subjects The number of subjects on each arm 
##' @param number.dropouts The number of subjects who dropout on each arm
##' @param total.events The total number of events for each arm
##' @param time.at.risk The total time at risk for each arm
##' @param empirical.rates total.events/time.at.risk 
##' 
##' The \code{print.summary.SingleSim} method has been implemented
##' 
##' @name summary.SingleSim.object
NULL
