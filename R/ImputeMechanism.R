##' @export
print.ImputeMechanism <- function(x,...){
  cat("Imputation Method:",x$name)
  .internal.output.list(x$parameters)
}


##' @export
GetImputedDataSet <- function(imputeSim,index){
  
  ValidateGetImputeDSArgs(imputeSim,index)
  
  retVal <- imputeSim$fit$singleSim
  retVal$status <- "imputed"
  retVal$data$censored.time <- imputeSim$imputed.values[,index]$new.censored.times
  retVal$event.times <- mapply(c,retVal$event.times,imputeSim$imputed.values[,index]$newevent.times,SIMPLIFY = FALSE)
  retVal$data$observed.events <-  vapply(retVal$event.times,length,FUN.VALUE=numeric(1))
  return(retVal)
}



.internal.number.data.sets <- function(imputedSim){
  return(ncol(imputedSim$imputed.values))
}

ValidateGetImputeDSArgs <- function(imputeSim,index){
  if(class(imputeSim)!="ImputeSim"){
    stop("Invalid argument: imputeSim argument must be an ImputeSim object")
  }
  
  if(!is.numeric(index) || is.na(index) || index < 0 || length(index)>1 || !.internal.is.wholenumber(index)){
    stop("Invalid argument: index")
  }
  
  if(index > .internal.number.data.sets(imputeSim)){
    stop("index too big, not enough imputed data sets!")
  }
}


.internal.impute <- function(fit,treatment.p.choice){
  
  df <- fit$singleSim$data
  
  lapply(1:nrow(df), function(i){
    study.time <- fit$singleSim$study.time 
    time.left <- study.time - df$censored.time[i]
    if(time.left==0){return(numeric(0))}
    
    if(df$arm[i]==0){
      p <- (fit$p[1] * time.left)/(fit$gamma[1]+fit$p[1]*study.time)
      gamma <- fit$gamma[1] + df$observed.events[i]
    }
    else{
      p <- (fit$p[1]*time.left)/(fit$gamma[2] + fit$p[2]*df$censored.time[i] + fit$p[1]*time.left)
      p <- c(p, (fit$p[2]*time.left)/(fit$gamma[2]+fit$p[2]*study.time))
      p <- treatment.p.choice(p)
      gamma <- fit$gamma[2] + df$observed.events[i]
    }  
      
    u <- (p/(1-p))*gamma
    
    rate <- GetSimRates(time.left,number.subject=1,event.rate=u/time.left,dispersion=1/gamma) 
    return(GetEventTimes(rate,time.left)+df$censored.time[i])  
    #note numeric(0)+x = numeric(0)
  })
  
}