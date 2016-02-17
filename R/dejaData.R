##' Boo
##' @name DejaData.object
NULL

##' Create a \code{DejaData} object
##' 
##' This object is can be used to create a \code{SingleSim} object with
##' subject specific rates
##' 
##' @param data A data frame containing the subject
##' @param arm, character the column name of the treatment arm for each subject
##' @param rate, character the column name of the rate to be used when simulating
##' @param Id, character the column name of subject Id
##' @export
MakeDejaData <- function(data,arm,rate,Id){
  
  #validateStuff
    
  retVal <- list(data=data,
       arm=arm,
       rate=rate,
       Id=Id)
  
  class(retVal) <- "DejaData"
  retVal
}


defaultDejaData <- function(number.subjects,event.rates){
  expandtoPerSubject <- function(val,number.subjects){
    unlist(mapply(rep,val,number.subjects,SIMPLIFY = FALSE))  
  }

  MakeDejaData(data=data.frame(arm=expandtoPerSubject(c(0,1),number.subjects),
                               rate=expandtoPerSubject(event.rates,number.subjects),
                               Id=1:sum(number.subjects)),
                               arm="arm",rate="rate",Id="Id")
}  
  