##' ImputeMechanism object
##' 
##' An object which defines a mechanism for taking
##' a \code{SingleSimFit} object and imputing missing data
##' to create a \code{ImputeSim}
##' 
##' It is possible to create user defined mechanisms, however, 
##' common mechanisms have already been implemented. For example see
##' \code{\link{weighted_j2r}} 
##' 
##' A \code{print.ImputeMechanism} method is defined.
##'
##' @section Structure: The following components must be included in
##' an ImputeMechanism Object
##' @param name The method name (used for printing)
##' @param cols.needed which columns of the SingleSim data frame are required by the method, typically
##' \code{c("censored.time","observed.events","arm")}
##' @param impute A function which takes a \code{SingleSimFit} object and outputs the details for a single
##' imputed data set, specifically a list with two elements:
##' \code{new.censored.times} - a vector of times subjects were censored (after taking into account imputation)
##' and \code{newevent.times} - a list of vectors where the vectors contain the imputed event times for the subjects
##' (these vectors do not contain the observed event times before subject drop out). If a subject has no imputed events then 
##' the vector \code{numeric(0)} is returned.   
##' @param parameters A list of named parameters describing the method (used for printing) - or NULL if none
##' @examples 
##' j2r <- weighted_j2r(trt.weight=0)
##' @name ImputeMechanism.object
NULL


##' @export
print.ImputeMechanism <- function(x,...){
  cat("Imputation Method:",x$name)
  .internal.output.list(x$parameters)
}

