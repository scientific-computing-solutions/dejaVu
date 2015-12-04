##' Scenario object
##' 
##' TODO
##' @name Scenario.object
NULL 


##' @export
as.data.frame.Scenario <- function(x,row.names = NULL, optional = FALSE,use.adjusted.pval=FALSE,...){
  
  Validate.adjusted.pval(x$summaries,use.adjusted.pval)
  
  .extract <- function(name,fun.val=numeric(1)){
    vapply(x$summaries,function(y){y[[name]]},FUN.VALUE = fun.val)
  }
  
  pval.string <-if(use.adjusted.pval)  "adjusted.pval" else "pval"
  
  
  data.frame(replica=1:length(x$summaries),
             treatment.effect=.extract("treatment.effect"),
             se=.extract("se"),
             pval=.extract(pval.string),
             theta=.extract("theta")
            )

}


Validate.adjusted.pval <- function(x,use.adjusted.pval){
  if(!is.logical(use.adjusted.pval) || length(use.adjusted.pval)>1){
   stop("Invalid argument use.adjusted.pval")
  }

  if(use.adjusted.pval && class(x[[1]])=="summary.SingleSimFit") {
    stop("Invalid arguments: cannot use adjusted.pval if multiple imputed test statistic is not used")
  }
}


##' Boo
##' @param object TODO
##' @param description TODO
##' @return TODO
##' @export
CreateScenario <- function(object,description=""){

  #validate description
  
  if(!is.list(object)){
    stop("Invalid argument: object must be a list of either summary.ImputeSimFit or summary.SingleSimFit")
  }    
  
  classes <- lapply(object,class)
  if(!all(classes=="summary.ImputeSimFit") && !(all(classes=="summary.SingleSimFit"))){
    stop("Invalid argument: object must be a list of either summary.ImputeSimFit or summary.SingleSimFit")
  }
  
  retVal <- list(description=description,
                 summaries=object)
  class(retVal) <- "Scenario"                             
  retVal                             
}


##' summary.Scenario object
##' 
##' TODO
##' @name summary.Scenario.object
NULL 

##' @export
summary.Scenario <- function(object,alpha=0.05,use.adjusted.pval=FALSE,...){
  Validate.adjusted.pval(object$summaries,use.adjusted.pval)
  
  data <- as.data.frame.Scenario(object,use.adjusted.pval)
  
  retVal <- list(treatment.effect=exp(mean(log(data$treatment.effect))),
                 se=mean(data$se),
                 power=sum(data$pval<alpha)/nrow(data),
                 alpha=alpha,
                 use.adjusted.pval=use.adjusted.pval,
                 description=object$description)
  
  class(retVal) <- "summary.Scenario"
  retVal 
}


##' @export
print.summary.Scenario <- function(x,...){
  cat("Scenario summary\n")
  if(!x$description==""){
    cat("Description:",x$description,fill=TRUE)
  }
  cat("Estimated Treatment Effect:",x$treatment.effect,fill=TRUE)
  cat("Event Rate Reduction: ",100*(1-x$treatment.effect) ,"%\n",sep="")
  cat("Estimated SE (log) treatment effect:",x$se,fill=TRUE)
  cat("Using alpha=",x$alpha,sep="")
  if(x$use.adjusted.pval){
    cat("(and an adjusted number of degrees of freedom)")
  }
  cat("\n  power:",x$power,fill=TRUE)
}