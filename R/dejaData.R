##' Data frame of covariates for simulating recurrent events
##' 
##' This object allows covariates to be included in the simulation procedure
##' The object is created using the \code{\link{MakeDejaData}} function
##' 
##' @section Structure: The above components must be included in
##' a DejaData Object
##' @param data A data frame containing the subject
##' @param arm, character the column name of the treatment arm for each subject
##' @param rate, character the column name of the rate to be used when simulating
##' @param Id, character the column name of subject Id 
##' @name DejaData.object
NULL

##' Create a \code{DejaData} object
##' 
##' This object is can be used to create a \code{SingleSim} object with
##' subject specific rates
##' 
##' @param data A data frame containing the subject
##' @param arm, character the column name of the treatment arm for each subject
##' @param Id, character the column name of subject Id
##' @param rate, character the column name of the rate to be used when simulating (or NULL,
##' if using DejaData to import a data set, see \code{\link{ImportSim}})
##' @return A \code{DejaData} object
##' @examples
##' 
##' set.seed(232)
##' 
##' my.df <- data.frame(Id=1:100,
##'                     arm=c(rep(0,50),rep(1,50)),
##'                     covar=rbinom(n=100,size=1,prob=0.5))
##' 
##' my.df$rate <- 0.0025 + my.df$covar*0.002 + (1-my.df$arm)*0.002
##' 
##' my.dejaData <- MakeDejaData(my.df,arm="arm",rate="rate",Id="Id")
##' 
##' @export
MakeDejaData <- function(data,arm,Id,rate=NULL){
  
  if(!"data.frame" %in% class(data)){
    stop("data must be a data frame")
  }
  
  if(any(!c(arm,Id) %in% colnames(data))){
    stop("Invalid argument check the arm, rate and Id column names")
  }
  
  if(nrow(data)<2){
    stop("Data set must contain at least two rows")
  }
    
  if(any(!data[,arm] %in% c(0,1))){
    stop("Treatment arm must be '0' or '1'")
  }
  
  if(any(!c(0,1) %in% data[,arm])){
    stop("Data frame only contains one treatment group!")
  }
  
  if(!is.null(rate)){
    if(!rate %in% colnames(data)){
      stop("Invalid rate column")
    }
    if(any(!is.numeric(data[,rate]) || data[,rate] < 0)){
      stop("Invalid rate, must be non-negative number")
    }
  }  
  
  if(nrow(data)!=length(unique(data[,Id]))){
    stop("Subject Ids must be unique")
  }
  
  retVal <- list(data=data,
       arm=arm,
       rate=rate,
       Id=Id)
  
  class(retVal) <- "DejaData"
  retVal
}


#function which creates a DejaData object for a given number of 
#subjects per arm and per arm event.rates, where no covariates
#are used
defaultDejaData <- function(number.subjects,event.rates){
 
  validateDefaultDejaDataArgs(number.subjects,event.rates)
  
  if(length(number.subjects)==1) number.subjects <- rep(number.subjects,2)
  
   expandtoPerSubject <- function(val,number.subjects){
    unlist(mapply(rep,val,number.subjects,SIMPLIFY = FALSE))  
  }

  MakeDejaData(data=data.frame(arm=expandtoPerSubject(c(0,1),number.subjects),
                               rate=expandtoPerSubject(event.rates,number.subjects),
                               Id=1:sum(number.subjects)),
                               arm="arm",rate="rate",Id="Id")
}  
  