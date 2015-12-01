ValidateSimCompleteArgs <- function(study.time,number.subjects,event.rates,dispersions){
  return(NULL)
}


GetSimRates <- function(study.time,number.subject,event.rate,dispersion){
  if(dispersion==0){
    return(rep(event.rate,number.subject))  
  }
  
  mean.event <- event.rate*study.time
  p <- (dispersion*mean.event)/(1+dispersion*mean.event)
  
  lambda <- rgamma(n=number.subject,shape=1/dispersion,scale=p/(1-p))
  return(lambda/study.time)
}



# does not take into account the requirement for the adjacent 
# exacerbations to be > 7 days apart
GetEventTimes <- function(rate,study.time){
  event.times <- numeric(0)
  current.time <- rexp(1,rate)
  
  while(current.time <= study.time){
    event.times <- c(event.times,current.time)
    current.time <- current.time + rexp(1,rate)
  }
  
  return(event.times)
}
