##' ImputeMechanism object
##' 
##' An object which defines a mechanism for taking
##' a \code{SingleSimFit} object and imputing missing data
##' to create a \code{ImputeSim}
##' 
##' It is possible to create user defined mechanisms, however, certain
##' common mechanisms have already been implemented. For example see
##' \code{\link{j2r}} 
##' 
##' A \code{print.ImputeMechanism} method is defined.
##'
##' @section Structure: The following components must be included in
##' an ImputeMechanism Object
##' @param name 
##' @param cols.needed which columns of the SingleSim data frame are required by the method, typically
##' \code{c("censored.time","observed.events","arm")}
##' @param impute A function which takes a \code{SingleSimFit} object and outputs the details for a single
##' imputed data set, specifically a list with two elements:
##' \code{new.censored.times} - a vector of times subjects were censored (after taking into account imputation)
##' and \code{newevent.times} - a list of vectors where the vectors contain the imputed event times for the subjects
##' (these vectors do not contain the observed event times before subject drop out). If a subject has no imputed events then 
##' the vector \code{numeric(0)} is returned.   
##' @param parameters A list of named parameters describing the method (used for printing) - or NULL if none
##' @examples 
##' j2r.procedure <- j2r()
##' @name ImputeMechanism.object
NULL


##' @export
print.ImputeMechanism <- function(x,...){
  cat("Imputation Method:",x$name)
  .internal.output.list(x$parameters)
}


.internal.impute <- function(fit,treatment.p.choice){
  
  df <- fit$singleSim$data
  
  lapply(1:nrow(df), function(i){
    study.time <- fit$singleSim$study.time 
    time.left <- study.time - df$censored.time[i]
    if(time.left==0){return(numeric(0))}
    
    if(df$arm[i]==0){
      p <- (fit$p[1] * time.left)/(fit$gamma[1]+fit$p[1]*study.time)
      gamma <- fit$gamma[1] + df$observed.events[i]
    }
    else{
      p <- (fit$p[1]*time.left)/(fit$gamma[2] + fit$p[2]*df$censored.time[i] + fit$p[1]*time.left)
      p <- c(p, (fit$p[2]*time.left)/(fit$gamma[2]+fit$p[2]*study.time))
      p <- treatment.p.choice(p)
      gamma <- fit$gamma[2] + df$observed.events[i]
    }  
      
    u <- (p/(1-p))*gamma
    
    rate <- GetSimRates(time.left,number.subject=1,event.rate=u/time.left,dispersion=1/gamma) 
    return(GetEventTimes(rate,time.left)+df$censored.time[i])  
    #note numeric(0)+x = numeric(0)
  })
  
}