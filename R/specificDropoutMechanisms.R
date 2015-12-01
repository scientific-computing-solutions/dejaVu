##' TODO
##' @export
ConstantRateDrop <- function(log.rate,var=0){
  
  #validation here
  
  text <- paste("Constant Dropout Rate","TODO") 
  cols.needed <- "censored.time"  
  
  sd <- sqrt(var)
  
  f <- function(event.times,data){
    rate <- exp(log.rate+rnorm(1,mean = 0,sd = sd))
    dropout.time <- rexp(1,rate)
    return(pmin(dropout.time,data$censored.time))
  }
    
  retVal <- list(
    type="MCAR",
    text=text,
    cols.needed=cols.needed,
    GetDropTime=f
  )
  
  class(retVal) <- "DropoutMechanism"
  return(retVal)
}