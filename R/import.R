##' Import an existing data frame for use with the package
##'
##' @param dejaData a \code{DejaData} object contain the subject
##' covariates and treatment arm
##' @param event.times A list of vectors, containing the observed
##' event times of each subject. If no events are observed then
##' numeric(0) should be used. See example in this help file for more
##' details
##' @param status The status of the data set imported, either
##' "complete" (if all subjects complete their follow up period) or
##' "dropout" (if not)
##' @param study.time The total follow up time according to study
##' protocol
##' @param censored.time If status is "dropout", this is a vector of
##' the times at which each subject is censored
##' @param actual.events If status is "dropout" and the total number
##' of events (i.e. not just the number observed) is known (e.g. if a
##' different simulation procedure was used) a vector of total number
##' of events should be included. If the number is not known or status
##' is "complete" then this should be set to NULL
##' @param allow.beyond.study Whether or not to allow imported data
##' with events after the nominal end of study.
##' @return A SingleSim object   
##' @examples 
##' 
##' covar.df <- data.frame(Id=1:6,
##'                        arm=c(rep(0,3),rep(1,3)),
##'                        Z=c(0,1,1,0,1,0))
##'  
##' dejaData <- MakeDejaData(covar.df,arm="arm",Id="Id") 
##' 
##' 
##' event.times <- list(c(25,100,121,200,225),
##'                     c(100,110),c(55),numeric(0),
##'                     150,45)     
##'
##' complete.dataset <- ImportSim(dejaData, event.times,
##'                     status="complete",
##'                     study.time=365)
##'  
##' censored.time  <- c(365,178,100,245,200,100)
##' 
##' dropout.dataset <- ImportSim(dejaData, event.times,
##'                     status="dropout",
##'                     study.time=365,
##'                     censored.time=censored.time)
##' 
##' @export
ImportSim <- function(dejaData, event.times, status, study.time, 
                      censored.time=NULL, actual.events=NULL,
                      allow.beyond.study=FALSE){
 
  if(!status %in% c("complete","dropout")){
    stop("Invalid status, it must be 'complete' or 'dropout'")
  }
  
  if(!.internal.is.finite.number(study.time) || study.time < 0){
    stop("Invalid study.time argument it must be a single positive finite numeric value")
  }
  
  if(class(dejaData)!= "DejaData"){
    stop("dejaData must be a DejaData object")
  }
  
  if(status=="complete" &&
     (!is.null(actual.events) ||
          !is.null(censored.time))){
      stop("Cannot give actual.events or censored.time arguments if status is complete")
  }
  
  if(status=="dropout" && is.null(censored.time)){
    stop("Must give censored.time if status is dropout")
  }
  
  number.subjects <- nrow(dejaData$data)
  importvalidate.event.times(event.times,
                             study.time,
                             number.subjects,
                             allow.beyond.study) 
  observed.events <-
      vapply(event.times,
             length,
             FUN.VALUE=numeric(1))
  
  
  if(status=="complete"){
      actual.events <- observed.events
      ## we have validated event times
      ## which are numbers
      ## so now let's make fake censored times
      censored.time <-
          vapply(event.times,
                 function(which) {max(study.time, which)},
                 0)
  }
  
  if(is.null(actual.events)){
    actual.events <- rep(as.numeric(NA),number.subjects)
  }  
  
  importvalidate.actual.events(actual.events, observed.events)        
  importvalidate.censored.time(censored.time, study.time, event.times,
                               allow.beyond.study)  
  
  data <- data.frame(Id=dejaData$data[,dejaData$Id],
                     arm=as.factor(dejaData$data[,dejaData$arm]),
                     censored.time=censored.time,
                     observed.events=observed.events,
                     actual.events=actual.events)
  
  .singleSimConstructor(data=cbind(data, remove.dejacols(dejaData)),
                        event.times=event.times,
                        status=status,
                        subject.rates=NULL,  
                        dropout.mechanism=if(status=="dropout") "User imported dropout data" else NULL,
                        impute.mechanism=NULL,
                        study.time=study.time,
                        event.rates=NULL,
                        dispersions=NULL)

}


#check that event.times are a list of the appropriate length and each element is a 
#sorted vector of numbers >= 0
importvalidate.event.times <- function(event.times, study.time,
                                       number.subjects,
                                       allow.beyond.study){
  if(!is.list(event.times) || length(event.times) != number.subjects){
    stop("event.times must be a list of length",number.subjects,"\n")
  }

  checker <-
      if (isTRUE(allow.beyond.study)) {
          function(x) {
              if(!is.numeric(x))
                  stop("event.times must be numeric")
              if(any(x<0))
                  stop("Cannot have event.times < 0")
              if(length(x)!= 0 && sort(x) != x)
                  stop("Event times for subjects must be given in ascending order")
          }
      } else {
            function(x) {
                if(!is.numeric(x))
                    stop("event.times must be numeric")
                if (any((x<0) | (x > study.time)))
                    stop("Cannot have event.times < 0 or > study.time")
                if(length(x)!= 0 && sort(x) != x)
                    stop("Event times for subjects must be given in ascending order")
            }
        }
  
  lapply(event.times, checker)
}

#Check that actual.events is a vector of the appropriate length and all elements are >=
#observed events (which has been calculated from a validated event.times)
importvalidate.actual.events <- function(actual.events, observed.events){
  if(length(actual.events)!= length(observed.events)){
    stop("invalid length of actual.events")
  }

  lapply(actual.events,function(x){
    if(!is.na(x) && !!.internal.is.finite.number(x) &&
       !.internal.is.wholenumber(x)){
      stop("actual.events must be integers")
    }
  })
  
  
  if(any(!is.na(actual.events) & actual.events < observed.events)){
    stop("actual.events cannot be < the number of observed events")
  }
  
} 

#check that the censored.time argument is of correct length <= study.time and
#all event.times events occur before the censored time, note event.times
#has already been validated
importvalidate.censored.time <- function(censored.time,
                                         study.time,
                                         event.times,
                                         allow.beyond.study){
  
  if(length(censored.time) != length(event.times)) {
      stop("censored.time does not have one element per subject")
  }

  checker <-
      if (isTRUE(allow.beyond.study)) {
          function(x){
              if(!.internal.is.finite.number(x) || 
                 x < 0 ){
                  stop("censored.time must be non-negative")
              }
          }
      } else {
            function(x) {
                if(!.internal.is.finite.number(x) || 
                   (x < 0) ||
                   (x > study.time)){
                    stop("censored.time must be non-negative and <= study.time")
                }
            }
        }
  
  lapply(censored.time, checker)
  
  checkmax <- function(censored.time,event.times){
    if (censored.time < max(-Inf, event.times)) {
      stop("event.times cannot occur after censored.time")
    }
  }
  
  mapply(checkmax, censored.time, event.times)
  
}  
