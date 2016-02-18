##' Import an existing data frame for use with the package
##'
##' @param dejaData a \code{DejaData} object contain the subject covariates and treatment arm
##' @param event.times A list of vectors, containing the observed event times of each subject. If no events 
##' are observed then numeric(0) should be used. See example in this help file for more details 
##' @param status The status of the data set imported, either "complete" (if all subjects complete their 
##' follow up period) or "dropout" (if not) 
##' @param study.time The total follow up time according to study protocol
##' @param censored.time If status is "dropout", this is a vector of the times at which each subject is censored  
##' @param actual.events If status is "dropout" and the total number of events (i.e. not just the number observed)
##' is known (e.g. if a different simulation procedure was used) a vector of total number of events should be included. If
##' the number is not known or status is "complete" then this should be set to NULL 
##' known (this is real data and a subject drops out) then use NA
##' @return A SingleSim object   
##' @examples 
##' 
##' set.seed(10)
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
                      censored.time=NULL, actual.events=NULL){
 
  if(!status %in% c("complete","dropout")){
    stop("Invalid status, it must be 'complete' or 'dropout'")
  }
  
  if(!.internal.is.finite.number(study.time) || study.time < 0){
    stop("Invalid study.time argument it must be a single positive finite numeric value")
  }
  
  if(class(dejaData)!= "DejaData"){
    stop("dejaData must be a DejaDate object")
  }
  
  if(status=="complete" && (!is.null(actual.events) || !is.null(censored.time))){
    stop("Cannot give actual.events or censored.time arguments if status is complete")
  }
  
  if(status=="dropout" && is.null(censored.time)){
    stop("Must give censored.time if status is dropout")
  }
  
  number.subjects <- nrow(dejaData$data)
  importvalidate.event.times(event.times,study.time,number.subjects) 
  observed.events <- vapply(event.times,length,FUN.VALUE=numeric(1))
  
  
  if(status=="complete"){
    actual.events <- observed.events
    censored.time <- rep(study.time, number.subjects)  
  }
  
  if(is.null(actual.events)){
    actual.events <- rep(NA,number.subjects)
  }  
  
  importvalidate.actual.events(actual.events,observed.events)        
  importvalidate.censored.time(censored.time,event.times)  
  
  data <- data.frame(Id=dejaData$data[,dejaData$Id],
                     arm=as.factor(dejaData$data[,dejaData$arm]),
                     censored.time=censored.time,
                     observed.events=observed.events,
                     actual.events=if(status=="complete") observed.events else rep(as.numeric(NA),number.subjects))
  
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


importvalidate.event.times <- function(event.times,study.time,number.subjects){
  
}

importvalidate.actual.events <- function(actual.events,observed.events){
  
} 

importvalidate.censored.time <- function(censored.time,event.times){
  
}  