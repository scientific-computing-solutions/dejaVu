##' DropoutMechanism object
##' 
##' An object which defines a specific mechanism which  
##' takes a complete \code{SingleSim} object and returns
##' a set of drop out times for subjects
##' 
##' It is possible to create user defined mechanisms, however, certain
##' common mechanisms have already been implemented. For example see
##' \code{\link{ConstantRateDrop}} and \code{\link{LinearRateChangeDrop}}
##' 
##' Only the GetDropTime and cols.needed entries are required for
##' calculation, the other entries are used for printing the object 
##' 
##' \code{print.DropoutMechanism} methods is defined.
##'
##' @section Structure: The following components must be included in
##' a DropoutMechanism Object
##' @param type The type of mechanism (e.g. "MCAR" or "MNAR")
##' @param text A short string describing the mechanism (only used for printing)
##' @param cols.needed Which columns in the SingleSim$data data frame must be included for this
##' drop out mechanism to work. This option could allow drop out mechanism which depend on covariates to be
##' included. 
##' @param GetDropTime A function with two arguments event.times and data, the corresponding entries from the 
##' SingleSim object. This function should return a list of dropout times (if a subject does not dropout its dropout time should be 
##' their current censored.time (i.e. the study follow up time)) 
##' @param parameters A list of named parameters for the mechanism (only used for printing) or NULL if none
##' @name DropoutMechanism.object
##' @aliases print.DropoutMechanism 
NULL

##' @export
print.DropoutMechanism <- function(x,...){
  cat("Dropout Type:",x$type,fill=TRUE)
  cat("Dropout Mechanism:",x$text,fill=TRUE)
  .internal.output.list(x$parameters)
}



##' A function which creates a DropOut Mechanism object
##'
##' @inheritParams DropoutMechanism.object 
##' @seealso \code{\link{DropoutMechanism.object}}
##' @export
CreateNewDropoutMechanism <- function(type,text,cols.needed=vector("character"),GetDropTime,parameters=NULL){
  
  if(length(type)>1 || !type %in% c("MCAR","MAR","MNAR")){
    stop("Invalid type argument")
  } 
  
  if(length(text)>1){
    stop("Invalid text argument")
  }
  
  if(!is.vector(cols.needed)){
    stop("Invalid cols.needed argument should be a vector, if no columns needed then leave out argument")
  }
  
  if(!is.null(parameters) && !is.list(parameters)){
    stop("Invalid parameter list, if no parameters are needed, leave this argument NULL")
  }
  
  if(!is.function(GetDropTime)){
    stop("Invalid GetDropTime, it should be a function")
  }
  
  retVal <- list(type=type,
                 text=text,
                 cols.needed=cols.needed,
                 GetDropTime=GetDropTime,
                 parameters=parameters)
  
  class(retVal) <- "DropoutMechanism"
  retVal
}
