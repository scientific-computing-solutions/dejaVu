##' Create a Dropout Mechanism with constant dropout rate
##' 
##' Creates an MCAR \code{DropoutMechanism} object where subject \code{i} dropout
##' is exponentially distributed with rate \code{Ri} where 
##' \code{log(Ri) = C + Xi} for constant \code{C} and 
##' \code{Xi} a random normal variable with mean 0 and standard deviation \code{sigma}
##' 
##' @param log.rate \code{C} described in the details
##' @param var \code{sigma^2} described in the details section, by default = 0
##' @return A \code{DropoutMechanism} object
##' @seealso \code{\link{DropoutMechanism.object}}
##' @examples
##' ConstantRateDrop(log.rate=-6)
##' ConstantRateDrop(log.rate=-6,var=1)
##' @export
ConstantRateDrop <- function(log.rate,var=0){
  
  if(!is.numeric(log.rate) || is.na(log.rate) || length(log.rate)>1){
    stop("Invalid argument log.rate")
  }
  if(!is.numeric(var) || is.na(var) || length(var)>1 || var < 0){
    stop("Invalid argument var")
  }
  
  
  text <- "Constant rate per subject" 
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
    GetDropTime=f,
    parameters=list(log.rate=log.rate,between.subject.var=var)
  )
  
  class(retVal) <- "DropoutMechanism"
  return(retVal)
}