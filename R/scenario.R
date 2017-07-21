##' Scenario object
##' 
##' This class contains a collection of model fit summaries and summarizing
##' this object will calculate overall summary statistics such as power/type I error
##' 
##' Functions \code{as.data.frame.Scenario} and \code{summary.Scenario} have 
##' been implemented
##' 
##' @param description A string containind a description of the scenario
##' @param summaries A list of either \code{summary.ImputeSimFit} or \code{summary.SingleSimFit} objects
##' 
##' @seealso \code{\link{CreateScenario}} 
##' @name Scenario.object
NULL 


##' Extract the results of running a scenario
##' 
##' This function is a wrapper around \code{\link{CreateScenario}}
##' See the user guide vignette for an example of using this function
##' 
##' @param answer A named list of lists 
##' @param name The name of the lists of answer which should be extracted and 
##' put together into a \code{sc} 
##' @param description The description parameter to be passed into the \code{CreateScenario} function
##' @return A \code{Scenario} object  
##' @seealso \code{\link{CreateScenario}}
##' @export
extract_results <- function(answer,name,description){
  
  if(!is.character(name) || length(name) != 1){
    stop("Invalid argument name")
  }
  
  #using validation inside CreateScenario for more validation
  CreateScenario(lapply(answer,
                  function(x){
                    if(is.null(x[[name]])){
                      stop("Invalid name")
                    }   
                    x[[name]]}), description=description)  
}


##' @export
as.data.frame.Scenario <- function(x,row.names = NULL, optional = FALSE,use.adjusted.pval=FALSE,...){
  
  Validate.adjusted.pval(x$summaries,use.adjusted.pval)
  
  .extract <- function(name,fun.val=numeric(1)){
    vapply(x$summaries,function(y){y[[name]]},FUN.VALUE = fun.val)
  }
  
  dropout <- .extract("dropout",numeric(2))
  
  adjusted.string <- if(use.adjusted.pval) "adjusted." else ""
  
  data.frame(replica=1:length(x$summaries),
             treatment.effect=.extract("treatment.effect"),
             se=.extract("se"),
             pval=.extract(paste(adjusted.string,"pval",sep="")),
             df=.extract(paste(adjusted.string,"df",sep="")),
             dispersion=if(length(x$summaries[[1]]$dispersion)>0) .extract("dispersion") else NA,
             dropout.control.rate=dropout[1,]/x$summaries[[1]]$number.subjects[1],
             dropout.active.rate=dropout[2,]/x$summaries[[1]]$number.subjects[2]
            )

}

# validate the logical value use.adjusted.pval, it
# should be logical and must be FALSE if x is a list 
# of summary.SingleSimFit objects
Validate.adjusted.pval <- function(x,use.adjusted.pval){
  if(!is.logical(use.adjusted.pval) || length(use.adjusted.pval)>1){
   stop("Invalid argument use.adjusted.pval")
  }

  if(use.adjusted.pval && class(x[[1]])=="summary.SingleSimFit") {
    stop("Invalid arguments: cannot use adjusted.pval if multiple imputed test statistic is not used")
  }
}


##' Create \code{Scenario} object from list of Fit Summaries
##' @param object Either a list of \code{summary.SingleSimFit}
##' or \code{summary.ImputeSimFit} objects
##' @param description A character string describing the scenario 
##' (used for printing)
##' @return A \code{Scenario} object
##' @seealso \code{\link{Scenario.object}}
##' @export
CreateScenario <- function(object,description=""){

  if(length(description)!= 1){
    stop("Invalid description argument")
  }
  
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
##' This object contains the overall summary statistics for a specific
##' scenario. It is envisioned that multiple scenarios are run and a set of
##' \code{summary.Scenario} objects are created and these can then be used for
##' plotting
##' 
##' A \code{print.summary.Scenario} function has been implemented
##' 
##' @param treatment.effect The exp(mean(log(individual treatment effects))),
##' @param se The mean standard error of the (log) treatment effect
##' @param power The proportion of simulations for which the p-value is \code{< alpha}
##' @param alpha The significance level used when calculating power, by default 0.05 use
##' \code{summary(object,alpha=x)} to use a different p value
##' @param use.adjusted.pval logical, default FALSE should the p values calculated using
##' Rubin's formula with the adjusted number of degrees of freedom be used. Use \code{summary(object,use.adjusted.pval=TRUE)},
##' to use the adjusted p values
##' @param description A string containing a description of the scenario
##' @param dropout A list of summary statistics regarding number of subject dropouts 
##' @name summary.Scenario.object
NULL 

##' @export
summary.Scenario <- function(object,alpha=0.05,use.adjusted.pval=FALSE,...){
  Validate.adjusted.pval(object$summaries,use.adjusted.pval)
  if(!.internal.is.finite.number(alpha)|| alpha <=0 || alpha >=1){
    stop("Invalid alpha")
  }
  
  data <- as.data.frame.Scenario(object,use.adjusted.pval)
  
  retVal <- .internal.summary.Scenario(data,alpha)
  retVal$description <- object$description
  retVal$use.adjusted.pval <- use.adjusted.pval
  class(retVal) <- "summary.Scenario"
  retVal 
}


.internal.summary.Scenario <- function(data,alpha){
  list(treatment.effect=exp(mean(log(data$treatment.effect))),
    se=mean(data$se),
    power=sum(data$pval<alpha)/nrow(data),
    alpha=alpha,
    dropout=list(summary(data$dropout.control.rate),
                 summary(data$dropout.active.rate)))
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
    cat(" (and an adjusted number of degrees of freedom)")
  }
  cat("\n  power:",x$power,fill=TRUE)
  
  cat("Dropout rate summary statistics\n")
  cat("Control arm:\n")
  print(x$dropout[[1]])
  cat("Active arm:\n")
  print(x$dropout[[2]])
}