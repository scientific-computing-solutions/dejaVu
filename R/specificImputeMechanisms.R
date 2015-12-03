##' Create a j2r \code{ImputeMechanism} object
##' 
##' Imputation using this mechanism will follow the jump to reference
##' (j2r) model whereby missing counts for subjects in both arms will 
##' be imputed according to the mean of the placebo arm conditioned
##' on the subject's observed number of events
##' 
##' @return An \code{ImputeMechanism} object
##' @seealso \code{\link{ImputeMechanism.object}}
##' @export
j2r <- function(){
  
  treatment.p.choice <- function(ps){
    return(ps[1])
  }
  
  
  f <- function(fit){
    return(list(newevent.times=.internal.impute(fit,treatment.p.choice),
                new.censored.times=rep(fit$singleSim$study.time,numberSubjects(fit))
    ))
    
  }
  
  retVal <- list(name="j2r",
                 cols.needed=c("censored.time","observed.events","arm"),
                 impute=f,
                 parameters=NULL)
  class(retVal) <- "ImputeMechanism"
  return(retVal)
}

##' TODO here
##' @export
efficacy <- function(){
  
  treatment.p.choice <- function(ps){
    return(ps[2])
  }
  
  
  f <- function(fit){
    return(list(newevent.times=.internal.impute(fit,treatment.p.choice),
                new.censored.times=rep(fit$singleSim$study.time,numberSubjects(fit))
    ))
    
  }
  
  retVal <- list(name="efficacy",
                 cols.needed=c("censored.time","observed.events","arm"),
                 impute=f,
                 parameters=NULL)
  class(retVal) <- "ImputeMechanism"
  return(retVal)
}

##' TODO here
##' @export
tr <- function(){
  
  treatment.p.choice <- function(ps){
    return(mean(ps))
  }
  
  
  f <- function(fit){
    return(list(newevent.times=.internal.impute(fit,treatment.p.choice),
                new.censored.times=rep(fit$singleSim$study.time,numberSubjects(fit))
    ))
    
  }
  
  retVal <- list(name="tr",
                 cols.needed=c("censored.time","observed.events","arm"),
                 impute=f,
                 parameters=NULL)
  class(retVal) <- "ImputeMechanism"
  return(retVal)
}