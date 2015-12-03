##' SingleSimFit object
##' 
##' TODO
##' 
##' @name SingleSimFit.object
NULL


##' summary.SingleSimFit
##' 
##' TODO
##' @name summary.SingleSimFit
NULL

##' @export
summary.SingleSimFit <- function(x,...){
  if(!x$equal.dispersion){
    stop("Cannot generate a summary if equal.dispersion is FALSE")
  }  
 
  model.summary <- summary(x$model)
  
  retVal <- list(model.summary=model.summary,
                 rate.ratio=exp(model.summary$coefficient[2,1]),
                 se.corrected=model.summary$coefficient[2,2]*exp(model.summary$coefficient[2,1]),
                 se.logrates=model.summary$coefficient[2,2], #se2
                 theta=model.summary$theta, #only if negative binomial
                 exp.intercept=exp(model.summary$coefficient[1,1])*c(1,exp(model.summary$coefficient[2,1])),         
                 pval=model.summary$coefficient[2,4],
                 datastatus=x$singleSim$status
                 )# exp.intercept this was nbd.mp and nbd.mt in old code 
  
  class(retVal) <- "summary.SingleSimFit"
  retVal
}

##' @export
print.summary.SingleSimFit <- function(x,...){
  cat("Summary for model fit for",x$datastatus,"data set",fill = TRUE)
  cat("Ratio Rate:",x$rate.ratio,fill=TRUE)
  cat("SE corrected:",x$se.corrected,fill=TRUE)
  cat("SE logrates:",x$se.logrates,fill=TRUE)
  cat("exp intercept:",x$exp.intercept,fill=TRUE)
  if(!is.null(x$theta)){
    cat("Negative binomial model parameter theta:",x$theta,fill=TRUE)
  }
  cat("p-value:",x$pval,fill=TRUE)
  cat("Note p-value is associated with this individual data set\n")
}


##' Produce imputed data sets
##' 
##' Given a \code{SingleSimFit} object (with impute.parameters not NULL)
##' and an imputation mechanism,
##' create a collection of imputed data sets 
##' 
##' @param fit A \code{SingleSimFit} object
##' @param impute.mechanism An \code{ImputeMechanism} object
##' @param N The number of data sets to impute
##' @return An \code{ImputeSim} object
##' @export
Impute <- function(fit,impute.mechanism,N){
  validateImputeArguments(fit,impute.mechanism,N)

  retVal <- list(fit=fit,
                 impute.mechanism=impute.mechanism,
                 imputed.values=replicate(n=N, impute.mechanism$impute(fit),simplify="list"))
  class(retVal) <- "ImputeSim"
  return(retVal)
}

##' @export
numberSubjects.SingleSimFit <- function(x){
  numberSubjects(x$singleSim)
}

# check the arguments to Impute function are valid
validateImputeArguments <- function(fit,impute.mechanism,N){
  
  if(class(fit)!="SingleSimFit"){
    stop("Invalid fit argument, must be of class SingleSimFit")
  }
  
  lapply(fit$impute.parameters,function(x){
    if(is.null(x)) stop("Cannot impute using this SingleSimFit object (a negative binomial model was not fit)")
  })
  
  if(class(impute.mechanism)!="ImputeMechanism"){
    stop("Invalid mpute.mechanism argument, must one of class ImputeMechanism")
  }
  
  if(!is.numeric(N)|| is.na(N) || is.infinite(N) || length(N)>1 || N <= 0 || !.internal.is.wholenumber(N)){
    stop("Invalid argument N, must be positive integer")
  }
  
  
  if(!all(impute.mechanism$cols.needed %in% colnames(fit$singleSim$data))){
    stop("This impute mechanism requires ",paste(impute.mechanism$cols.needed,collapse=", "),
         "as column names in the simulated data frame")
  }
}
