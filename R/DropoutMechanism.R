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
##' @param parameters A list of named parameters for the mechanism (only used for printing)
##' @name DropoutMechanism.object
##' @aliases print.DropoutMechanism 
NULL

ValidColumns <- function(drop.mechanism,columnnames){
  if(!all(drop.mechanism$cols.needed %in% columnnames)){
    stop("This dropout mechanism requires ",paste(drop.mechanism$cols.needed,collapse=", "),
               "as column names in the simulated data frame")
  }
}

##' @export
print.DropoutMechanism <- function(x,...){
  cat("Dropout Type:",x$type,fill=TRUE)
  cat("Dropout Mechanism:",x$text,fill=TRUE)
  cat("Parameters:",fill=TRUE)
  invisible(lapply(seq_along(x$parameters),function(i){
    cat("  ",names(x$parameters)[i],": ",x$parameters[[i]],"\n",sep="")}))
}