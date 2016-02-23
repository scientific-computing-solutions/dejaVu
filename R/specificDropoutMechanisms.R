##' Create a Dropout Mechanism with constant dropout rate
##' 
##' Creates an MCAR \code{DropoutMechanism} object where subject \code{i} dropout
##' is exponentially distributed with rate \code{Ri} where 
##' \code{Ri = C*exp(Xi)} for constant \code{C} and 
##' \code{Xi} a random normal variable with mean 0 and standard deviation \code{sigma}
##' 
##' @param rate \code{C} described in the details
##' @param var \code{sigma^2} described in the details section, by default = 0
##' @return A \code{DropoutMechanism} object
##' @seealso \code{\link{DropoutMechanism.object}}
##' @examples
##' ConstantRateDrop(rate=0.0025)
##' ConstantRateDrop(rate=0.0025,var=1)
##' @export
ConstantRateDrop <- function(rate,var=0){
  
  if(!.internal.is.finite.number(rate) || rate <= 0){
    stop("Invalid argument rate")
  }
  if(!.internal.is.finite.number(var) || var < 0){
    stop("Invalid argument var")
  }
  
  
  text <- "Constant rate per subject" 
  cols.needed <- "censored.time"  
  sd <- sqrt(var)
  
  f <- function(event.times,data){
    rate <- rate*exp(rnorm(1,mean = 0,sd = sd))
    dropout.time <- rexp(1,rate)
    return(pmin(dropout.time,data$censored.time))
  }
  
  CreateNewDropoutMechanism(type="MCAR",text=text,
                            cols.needed=cols.needed,
                            GetDropTime=f,
                            parameters=list(rate=rate,between.subject.var=var)) 
  
}


##' Create a Dropout Mechanism with drop out rate which changes by a fixed constant
##' after every event
##' 
##' Creates an MAR \code{DropoutMechanism} object where subject \code{i} has piecewise exponential
##' dropout rate where the rate changes by a constant amount after each event, specifically after \code{j} events
##' the subject has rate \code{Rij = Cj*exp(Xij)} where \code{Cj=C+j*D} for constants \code{C}, \code{D} 
##' and Xij is a standard normal variable with mean 0 and standard deviation \code{sigma}   
##' 
##'  
##' @param starting.rate \code{C}, see description section.
##' @param rate.change \code{D}, see description section. Note if \code{D<0}, \code{Cj} could be negative
##' for large \code{j}, this is not possible and the rate remains constant if the next change would set \code{Cj<=0} 
##' @param var \code{sigma^2}, see description section
##' @return A \code{DropoutMechanism} object
##' @seealso \code{\link{DropoutMechanism.object}}
##' @examples
##' LinearRateChangeDrop(starting.rate=0.0025,rate.change=0.0005)
##' LinearRateChangeDrop(starting.rate=0.0025,rate.change=-0.00001,var=1)
##' @export
LinearRateChangeDrop <- function(starting.rate,rate.change,var=0){
  
  if(!.internal.is.finite.number(starting.rate) || starting.rate <= 0){
    stop("Invalid argument starting.rate")
  }
  if(!.internal.is.finite.number(var)|| var < 0){
    stop("Invalid argument var")
  }
  if(!.internal.is.finite.number(rate.change)){
    stop("Invalid argument rate.change")
  }
  
  
  text <- "Linear change in rate after each event" 
  cols.needed <- "censored.time"  
  sd <- sqrt(var)
  
  f <- function(event.times,data){
    
    current.time <- 0
    current.rate <- starting.rate
    event.times <- c(event.times,data$censored.time)
    end.index <- 1
    
    while(end.index <= length(event.times)){
      possible.drop <- current.time + rexp(1,current.rate*exp(rnorm(1,0,sd)))
      if(possible.drop < event.times[end.index]){
        return(possible.drop)
      }
      current.time <- event.times[end.index]
      current.rate <- current.rate + rate.change
      if(current.rate <= 0){
        current.rate <- current.rate - rate.change 
      }
      end.index <- end.index + 1
    }
    
    return(data$censored.time)  
    
  }
  
  CreateNewDropoutMechanism(type="MAR",text=text,
                            cols.needed=cols.needed,
                            GetDropTime=f,
                            parameters=list(starting.rate=starting.rate,
                                            rate.change.after.event=rate.change,
                                            between.subject.var=var)) 

}
