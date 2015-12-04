##' ImputeSim object
##' 
##' TODO
##' @name ImputeSim.object
NULL

##' TODO
##' 
##' @param imputeSim TODO
##' @param index TODO
##' @return TODO
##' @export
GetImputedDataSet <- function(imputeSim,index){
  
  ValidateGetImputeDSArgs(imputeSim,index)
  
  retVal <- imputeSim$fit$singleSim
  retVal$status <- "imputed"
  retVal$data$actual.censored.time <- retVal$data$censored.time 
  retVal$data$censored.time <- imputeSim$imputed.values[,index]$new.censored.times
  retVal$event.times <- mapply(c,retVal$event.times,imputeSim$imputed.values[,index]$newevent.times,SIMPLIFY = FALSE)
  retVal$data$observed.events <-  vapply(retVal$event.times,length,FUN.VALUE=numeric(1))
  retVal$impute.mechanism <- imputeSim$impute.mechanism
  return(retVal)
}



.internal.number.data.sets <- function(imputedSim){
  return(ncol(imputedSim$imputed.values))
}

ValidateGetImputeDSArgs <- function(imputeSim,index){
  if(class(imputeSim)!="ImputeSim"){
    stop("Invalid argument: imputeSim argument must be an ImputeSim object")
  }
  
  if(!.internal.is.finite.number(index) || index <= 0 || !.internal.is.wholenumber(index)){
    stop("Invalid argument: index")
  }
  
  if(index > .internal.number.data.sets(imputeSim)){
    stop("index too big, not enough imputed data sets!")
  }
}

##' @export
Simfit.ImputeSim <- function(x,family="negbin",equal.dispersion=TRUE,formula=GetDefaultFormula(equal.dispersion),...){
  if(!equal.dispersion){
    stop("Invalid argument equal.dispersion must be TRUE")
  }
  
  imputed.summaries <- lapply(1:.internal.number.data.sets(x),
                         function(index){
                           singleSim <- GetImputedDataSet(x,index)                     
                           summary(Simfit(singleSim,family=family,equal.dispersion=TRUE,formula=formula,...))  
                       })
  
  retVal <- list(imputeSim=x,
                 summaries=imputed.summaries)
  
  class(retVal) <- "ImputeSimFit"
  return(retVal)
}

##' @export
numberSubjects.ImputeSim <- function(x,...){
  numberSubjects(x$fit)
} 
