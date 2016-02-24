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
  
    
  dispersion <- if(length(x$summaries[[1]]$dispersion)!=0) .extract("dispersion") else NA
  
  data.frame(imputeID=1:.internal.number.data.sets(x$imputeSim),
             control.rate=rate.estimates[1,],
             active.rate=rate.estimates[2,],
             treatment.effect=.extract("treatment.effect"),
             se=.extract("se"),
             pval=.extract("pval"),
             dispersion=dispersion)  

}

  
##' @export
subjectsPerArm.ImputeSimFit <- function(x,...){
  subjectsPerArm(x$imputeSim)
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
##' @param number.subjects The number of subjects (per arm) for this imputed data set
##' @name summary.ImputeSimFit.object
NULL

##' @export
summary.ImputeSimFit <- function(object,...){
  
  data <- as.data.frame(object)

  N <- .internal.number.data.sets(object$imputeSim)
  
  if(N==1){
    stop("Cannot apply Rubin's formula for test statistic if only 1 imputed data set!")
  }
   
  retVal <- .rubinsformula(data$treatment.effect,data$se, object$summaries[[1]]$df,N)
  retVal$dropout <- object$imputeSim$dropout
  retVal$dispersion <- mean(data$dispersion)
  retVal$number.subjects <- subjectsPerArm(object)
  class(retVal) <- "summary.ImputeSimFit"
  retVal
}

.rubinsformula <- function(treatment.effects,ses,original.df,N){
  log.treatment.effects <- log(treatment.effects)
  se.log.treatment.effects <- ses
  
  Q <- mean(log.treatment.effects) #test statistic
  B <- var(log.treatment.effects) #between imputation var
  U <- mean(se.log.treatment.effects^2) #average imputation var
  
  se <- sqrt(U+(1+1/N)*B) #standard error of combined test statistic
  df <- (N-1)*(1+U/(B*(1+1/N)))^2 #(unadjusted d.o.f)
  
  v.0 <-original.df 
  v.hat <- (v.0*(v.0+1)/(v.0+3))*(1-(1+1/N)*B/(se*se))  
  adjusted.df <- 1/(1/df + 1/v.hat)
  
  getpval <- function(df){
    2*(1-pt(abs(Q/se),df=df))
  } 
  
  list(treatment.effect=exp(Q),
       se=se,
       df=df,
       adjusted.df=adjusted.df,
       pval=getpval(df=df),
       adjusted.pval=getpval(df=adjusted.df))
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