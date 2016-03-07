##' ImputeSim object
##' 
##' This object contains a collection of imputed data sets
##' derived from a \code{SingleSimFit} object and \code{ImputeMechanism}
##' 
##' @param singleSim The \code{SingleSim} object from which the imputed data sets
##' have been derived
##' @param impute.parameters The \code{SingleSimFit$impute.parameters} used to perform the imputation
##' @param impute.mechanism The \code{ImputeMechanism} object used to perform the
##' imputation
##' @param imputed.values A matrix with 1 column per imputed data set and two rows:
##' newevent.times a list of vectors containing the imputed event times (not including the events 
##' which were observed) and new.censored.times - a vector containing the times at which subjects (with imputed
##' data) are now censored
##' @param dropout A vector containing the number of subjects who have dropped out in each arm, for whom data is to be
##' imputed  
##' 
##' Use \code{\link{GetImputedDataSet}} to extract a single imputed data set and use \code{Simfit} to fit
##' a model to the set of data sets
##' 
##' @seealso \code{\link{GetImputedDataSet}}
##' @name ImputeSim.object
NULL

##' Output a single imputed data set
##' 
##' @param imputeSim A \code{ImputeSim} object which contains 
##' multiple imputed data sets
##' @param index numberic, which of the multiple imputed data sets to output
##' @return A \code{SingleSim} object with \code{status="imputed"}
##' @seealso \code{\link{ImputeSim.object}}
##' @export
GetImputedDataSet <- function(imputeSim,index){
  
  ValidateGetImputeDSArgs(imputeSim,index)
  
  retVal <- imputeSim$singleSim
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

# validation of GetImputedDatSet function
# see that function for parameter details
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
Simfit.ImputeSim <- function(x,family="negbin",equal.dispersion=TRUE,covar=NULL,...){
  if(!is.logical(equal.dispersion) || length(equal.dispersion)!=1 || !equal.dispersion){
    stop("Invalid argument equal.dispersion must be TRUE")
  }
  
  imputed.summaries <- lapply(seq_len(.internal.number.data.sets(x)),
                         function(index){
                           singleSim <- GetImputedDataSet(x,index)                     
                           summary(Simfit(singleSim,family=family,equal.dispersion=equal.dispersion,covar=covar,...))
                       })
  
  retVal <- list(imputeSim=x,
                 summaries=imputed.summaries)
  
  class(retVal) <- "ImputeSimFit"
  return(retVal)
}

##' @export
subjectsPerArm.ImputeSim <- function(x,...){
  subjectsPerArm(x$singleSim)
} 
