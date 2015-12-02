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
SimulateComplete <- function(study.time,number.subjects,event.rates,dispersions){
  
  ValidateSimCompleteArgs(study.time,number.subjects,event.rates,dispersions)
  
  if(length(number.subjects==1)){
    number.subjects <- rep(number.subjects,2)
  }
  
    
  subject.rates <- unlist(mapply(GetSimRates,study.time=study.time,
                          number.subject=number.subjects,
                          event.rate=event.rates,
                          dispersion=dispersions,SIMPLIFY=FALSE))

  event.times <- lapply(subject.rates,GetEventTimes,study.time=study.time)
  
  total.subjects <- sum(number.subjects)
  events <- vapply(event.times,length,FUN.VALUE=numeric(1))
  
  data <- data.frame(Id=1:total.subjects,
                     arm=as.factor(unlist(lapply(seq_along(number.subjects),function(x){rep(x-1,number.subjects[x])}))),   
                     censored.time=rep(study.time,total.subjects),
                     observed.events=events,
                     actual.events=events)
  
  retVal <- list(data=data,
                 event.times=event.times,
                 status="complete",
                 subject.rates=subject.rates, #needed? 
                 dropout.mechanism=NULL,
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
  }
  cat("Data:\n")
  str(x$data)
}


##' SingleSim Object
##' 
##' Stuff 
##'   
##' \code{print.SingleSim} and \code{summary.SingleSim} methods are defined.
##'
##' @section Structure: The following components must be included in
##' a SingleSim Object
##' @param data TODO
##' @param event.times TODO
##' @param status TODO
##' @param subject.rates TODO
##' @param dropout.mechanism TODO
##' @param study.time TODO
##' @param event.rates TODO
##' @param dispersions  TODO
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


##' S3 generic TODO
##' @param x TODO
##' @export
numberSubjects <- function(x){
  UseMethod("numberSubjects")
}

##' @export
numberSubjects.default <- function(x){
  stop("Invalid x for numberSubjects")
}

##' @export
numberSubjects.SingleSim <- function(x){
  return(nrow(x$data))  
}
