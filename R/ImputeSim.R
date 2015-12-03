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
  
  if(!is.numeric(index) || is.na(index) || index < 0 || length(index)>1 || !.internal.is.wholenumber(index)){
    stop("Invalid argument: index")
  }
  
  if(index > .internal.number.data.sets(imputeSim)){
    stop("index too big, not enough imputed data sets!")
  }
}
