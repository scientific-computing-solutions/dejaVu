##' ImputeSimFit object
##' 
##' An object which contains both a set of imputed data sets (\code{ImputeSim} object) and
##' a set of models fitted to them 
##' 
##' Calling \code{summary.ImputeSimFit} will apply Rubin's formula to calculate 
##' estimates for the treatment effect and standard error
##' 
##' Functions \code{summary.ImputeSimFit} and \code{as.data.frame.ImputeSimFit}
##' have been implemented
##' 
##' @param imputeSim The \code{ImputeSim} object for which models have been fitted
##' @param summaries A list of \code{summary.SingleSimFit} objects containing the model fits
##' for each of the imputed data sets
##' @seealso  \code{summary.ImputeSimFit} \code{\link{summary.SingleSimFit}} 
##' @name ImputeSimFit.object
NULL

##' @export
as.data.frame.ImputeSimFit <- function(x,row.names = NULL, optional = FALSE,...){
  
  .extract <- function(name,fun.val=numeric(1)){
    vapply(x$summaries,function(x){x[[name]]},FUN.VALUE = fun.val)
  }
  
  rate.estimates <- .extract("rate.estimate",fun.val=numeric(2))
  
    
  dispersion <- if(!is.null(x$summaries[[1]]$dispersion)) .extract("dispersion") else NA
  
  data.frame(imputeID=1:.internal.number.data.sets(x$imputeSim),
             control.rate=rate.estimates[1,],
             active.rate=rate.estimates[2,],
             treatment.effect=.extract("treatment.effect"),
             se=.extract("se"),
             pval=.extract("pval"),
             dispersion=dispersion)  

}

  
##' @export
numberSubjects.ImputeSimFit <- function(x,...){
  numberSubjects(x$imputeSim)
} 

##' summary.ImputeSimFit object
##' 
##' The summary of a \code{ImputeSimFit} object. Rubin's formula 
##' is used to combine the test statistics into a single summary
##' 
##' A \code{print.summary.ImputeSimFit} object has been implemented
##' 
##' @param treatment.effect The mean of the estimated treatment.effect
##' from the imputed data
##' @param se The standard error of the (log) treatment effect calculated using Rubin's formula
##' @param df The number of degrees of freedom used to calculate the p-value
##' @param adjusted.df The number of degrees of freedom used to calculate the adjusted
##' p-value (this should be used if the complete data number of degrees of freedom is small)
##' @param dispersion The mean of the estimated dispersion parameter 
##' @param pval The p-value for the test log(treatment.effect)=0 using Rubin's formula
##' @param adjusted.pval The p-value for the test log(treatment.effect)=0 using Rubin's
##' formula and the adjusted number of degrees of freedom 
##' @param dropout The number of subjects who drop out (per arm) for this imputed data set
##' @name summary.ImputeSimFit.object
NULL

##' @export
summary.ImputeSimFit <- function(object,...){
  
  data <- as.data.frame(object)
  
  #Rubin's formula 
  N <- .internal.number.data.sets(object$imputeSim)
  
  log.treatment.effects <- log(data$treatment.effect)
  se.log.treatment.effects <- data$se
  
  se <- sqrt((1+1/N)*var(log.treatment.effects)+ mean(se.log.treatment.effects^2))
  df <- (N-1)*(1+(mean(se.log.treatment.effects^2))/( (1+1/N)*var(log.treatment.effects)))^2
  
  v.0 <- object$summaries[[1]]$df
  v.hat <- (v.0*(v.0+1)/(v.0+3))*(1-(1+1/N)*(var(log.treatment.effects)/se^2))  
  adjusted.df <- 1/(1/df + 1/v.hat)
  
  getpval <- function(df){
    2*(1-pt(abs(mean(log.treatment.effects)/se) ,df=df))
  } 
 
  retVal <- list(treatment.effect=mean(data$treatment.effect),
                 se=se,
                 df=df,
                 adjusted.df=adjusted.df,
                 dispersion=mean(data$dispersion),
                 pval=getpval(df=df),
                 adjusted.pval=getpval(df=adjusted.df),
                 dropout=object$imputeSim$dropout
                )
  class(retVal) <- "summary.ImputeSimFit"
  retVal  
}

##' @export
print.summary.ImputeSimFit <- function(x,...){
  cat("Summary for imputed data sets",fill=TRUE)
  cat("Treatment Effect:",x$treatment.effect,fill=TRUE)
  cat("SE (log) treatment effect:",x$se,fill=TRUE)
  cat("Degrees of freedom:",x$df,fill=TRUE)
  cat("p-value:",x$pval,fill=TRUE)
  cat("Adjusted d.o.f:",x$adjusted.df,fill=TRUE)
  cat("Adjusted p-value:",x$adjusted.pval,fill=TRUE)
  if(!is.null(x$dispersion)){
    cat("Average dispersion:",x$dispersion,fill=TRUE)
  }
  cat("Number of subjects dropped out per arm:",x$dropout,fill=TRUE)
}