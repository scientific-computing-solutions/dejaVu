##' ImputeSimFit object
##' 
##' TODO
##' @name ImputeSimFit.object
NULL

##' @export
as.data.frame.ImputeSimFit <- function(x,...){
  
  .extract <- function(name,fun.val=numeric(1)){
    vapply(x$summaries,function(x){x[[name]]},FUN.VALUE = fun.val)
  }
  
  rate.estimates <- .extract("rate.estimate",fun.val=numeric(2))
  
  theta <- if(!is.null(x$summaries[[1]]$theta)) .extract("theta") else NULL
  
  data.frame(imputeID=1:.internal.number.data.sets(x$imputeSim),
             control.rate=rate.estimates[1,],
             active.rate=rate.estimates[2,],
             treatment.effect=.extract("treatment.effect"),
             se=.extract("se"),
             pval=.extract("pval"),
             theta=theta)  

}

  
##' @export
numberSubjects.ImputeSimFit <- function(x,...){
  numberSubjects(x$imputeSim)
} 

##' @export
summary.ImputeSimFit <- function(object,...){
  
  data <- as.data.frame(object)
  
  #Rubin's formula 
  N <- .internal.number.data.sets(object$imputeSim)
  
  log.treatment.effects <- log(data$treatment.effect)
  se.log.treatment.effects <- data$se
  theta <- data$theta

  se <- sqrt((1+1/N)*var(log.treatment.effects)+ mean(se.log.treatment.effects^2))
  df <- (N-1)*(1+(mean(se.log.treatment.effects^2))/( (1+1/N)*var(log.treatment.effects)))^2
  
  v.0 <- 2*numberSubjects(object)-2
  v.hat <- (v.0*(v.0+1)/(v.0+3))*(1-(1+1/N)*(var(log.treatment.effects)/se^2))  
  adjusted.df <- 1/(1/df + 1/v.hat)
  
  getpval <- function(df){
    2*(1-pt(abs(mean(log.treatment.effects)/se) ,df=df))
  } 
 
  retVal <- list(treatment.effect=mean(data$treatment.effect),
                 se=se,
                 df=df,
                 adjusted.df=adjusted.df,
                 theta=if(is.null(theta)) NULL else mean(theta),
                 pval=getpval(df=df),
                 adjusted.pval=getpval(df=adjusted.df)
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
  if(!is.null(x$theta)){
    cat("Average theta:",x$theta,fill=TRUE)
  }
  
}