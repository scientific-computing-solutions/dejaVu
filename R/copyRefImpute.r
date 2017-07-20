##' Create a copy reference \code{ImputeMechanism} object
##' 
##' Missing counts for subjects in both arms are imputed by assuming
##' the rate before and dropout are both equal to the control (reference) estimated rate.
##' This corresponds to what is usually termed the copy reference assumption.
##' 
##' @return An \code{ImputeMechanism} object
##' @seealso \code{\link{ImputeMechanism.object}}
##' @export
copy_reference <- function() {

  #A function which takes SimFit object and outputs a list of 2 elements
  #1) newevent.times a list of vectors, the imputed event times for each subject (if subject has no new imputed
  #events then the vector should be numeric(0))
  #2) new.censored.times - the time at which subjects are censored in the imputed data set
  f <- function(fit){
    return(list(newevent.times=.copy_ref_impute(fit),
                new.censored.times=pmax(fit$singleSim$data$censored.time,fit$singleSim$study.time)
    ))
    
  }
  
  CreateNewImputeMechanism(name="copy_reference",
                           cols.needed=c("censored.time","observed.events","arm"),
                           impute=f)
}


.copy_ref_impute <- function(fit){
  # performs the copy reference method using the given SimFit object
  # returns the imputed event times for each subjects as a list of vectors
  
  #@param fit is a SimFit object
  #@return returns a list of vectors, the imputed event times for each subject (if subject has no new imputed
  #events then the vector should be numeric(0))
  
  gamma_mu <- fit$genCoeff.function()
  
  #the data frame from the SimFit object  
  df <- fit$singleSim$data
  
  #take each subject in turn
  lapply(seq_len(nrow(df)), function(i){
    
    study.time <- fit$singleSim$study.time 
    time.left <- study.time - df$censored.time[i]
    if(time.left<=0){return(numeric(0))} #subject was not censored, so don't impute
    
    gamma <- gamma_mu$gamma[1] + df$observed.events[i]
    condrate <- ((gamma_mu$gamma[1]+df$observed.events[i])/
                   (gamma_mu$gamma[1]+gamma_mu$mu[i,1]*df$censored.time[i]))*gamma_mu$mu[i,1]
    #sanity check on this: when gamma parameter goes to zero (very high within-patient correlation),
    #the conditional rate reduces to the number of pre-dropout events observed for the subject
    #divided by their pre-dropout time.
    #When gamma goes to infinity, so their pre-dropout count has not predictive information
    #for their post dropout count, conditional rate reduces to marginal control rate
    
    rate <- GetSimRates(time.left,event.rate=condrate,dispersion=1/gamma) 
    return(GetEventTimes(rate,time.left)+df$censored.time[i])  
    #note numeric(0)+x = numeric(0)
  })
  
}