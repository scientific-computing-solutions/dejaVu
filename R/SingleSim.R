##' Simulate a complete data set
##' 
##' Simulate a complete data set of a recurrent event
##' clinical trial without dropouts using a negative binomial model
##' with given rates and dispersion parameters
##' 
##' TODO add some maths here?
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
                 subject.rates=subject.rates) #needed? 
  
  class(retVal) <- "SingleSim"
  return(retVal)
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
##' @name SingleSim.object
##' @aliases print.SingleSim summary.SingleSim
NULL


##' Simulate subject dropout
##' 
##' STUFF
##' 
##' @param simComplete TODO
##' @param drop.mechanism TODO
##' @return A \code{SingleSim} object with \code{status='dropout'}
##' @export 
SimulateDropout <- function(simComplete,drop.mechanism){
  validateSimulateDropout(simComplete,drop.mechanism)
  
  censored.time <- vapply(seq_along(simComplete$event.times),function(i){
    drop.mechanism$GetDropTime(event.times=simComplete$event.times[[i]],data=simComplete$data[i,])
  },FUN.VALUE=numeric(1))
  
  event.times <- mapply(function(censor.time,event.times){ 
                          return(event.times[event.times<=censor.time]) },
                        censor.time=censored.time,event.times=simComplete$event.times,SIMPLIFY = FALSE)
  
  simComplete$observed.events <- vapply(event.times,length,FUN.VALUE=numeric(1))
  simComplete$data$censored.time <- censored.time
  simComplete$status <- "dropout"
  return(simComplete)
}