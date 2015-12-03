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