##' @importFrom MASS glm.nb
NULL

##' S3 generic for fitting an negative binomial model
##' @param x The S3 object
##' @param equal.dispersion logical, should the arms have the same dispersion parameter when
##' fitting the negative binomial model
##' @param formula The formula to be used when calling \code{glm.nb} this should be
##' left as the default value for all but advanced users 
##' @param ... Additional arguments to be passed to the Surv function
##' @export
NBfit <- function(x,equal.dispersion,formula,...){
  UseMethod("NBfit")
}

##' @export
NBfit.default <- function(x,equal.dispersion,formula,...){
  stop("Invalid x for NBfit")
}


##' @export
NBfit.SingleSim <- function(x,equal.dispersion,formula=GetDefaultFormula(equal.dispersion)){
  
  data <- x$data[x$data$censored.time>0,]
  
  if(equal.dispersion){
    model <- glm.nb(formula=formula,data=data)
    mod.summary <- summary(model)
    gamma <- rep(mod.summary$theta,2)
    p <- exp(mod.summary$coeffi[1,1])*c(1,exp(mod.summary$coeffi[2,1]))
  }
  else{
    model <- lapply(0:1,function(i){glm.nb(formula=formula,data=data[data$arm==i,])})
    gamma <- vapply(model,function(mod){summary(mod)$theta},FUN.VALUE = numeric(1)) 
    p <-     vapply(model,function(mod){exp(summary(mod)$coeffi[1,1])},FUN.VALUE = numeric(1))
  }
  
  retVal <- list(singleSim=x,
                 equal.dispersion=equal.dispersion,
                 model=model,
                 gamma=gamma,
                 p=p)
  class(retVal) <- "SingleSimFit"
  return(retVal)
}


##' SingleSimFit object
##' 
##' TODO
##' 
##' @name SingleSimFit.object
NULL


##' @export
Impute <- function(fit,impute.mechanism,N){
  validateImputeArguments(fit,impute.mechanism,N)

  retVal <- list(fit=fit,
                 impute.mechanism=impute.mechanism,
                 imputed.values=replicate(n=N, impute.mechanism$impute(fit),simplify="list"))
  class(retVal) <- "ImputeSim"
  return(retVal)
}

numberSubjects.SingleSimFit <- function(x){
  numberSubjects(x$singleSim)
}
