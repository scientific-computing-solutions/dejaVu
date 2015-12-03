##' SingleSimFit object
##' 
##' TODO
##' 
##' @name SingleSimFit.object
NULL


##' Produce imputed data sets
##' 
##' Given a \code{SingleSimFit} object and an imputation mechanism,
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
