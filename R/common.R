##' @importFrom MASS glm.nb
NULL

#cat a named list
.internal.output.list <- function(x){ 
  if(!is.null(x)){
    cat("Parameters:",fill=TRUE)
    invisible(lapply(seq_along(x),function(i){
      cat("  ",names(x)[i],": ",x[[i]],"\n",sep="")}))
  }
}

.internal.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

##' S3 generic to output the number of subjects in a given object
##' @param x The object
##' @return The number of subjects
##' @export
numberSubjects <- function(x){
  UseMethod("numberSubjects")
}

##' @export
numberSubjects.default <- function(x){
  stop("Invalid x for numberSubjects")
}

##' S3 generic for fitting an negative binomial model
##' @param x The S3 object
##' @param equal.dispersion logical, should the arms have the same dispersion parameter when
##' fitting the negative binomial model
##' @param formula The formula to be used when calling \code{glm.nb} this should be
##' left as the default value for all but advanced users 
##' @param ... Additional arguments to be passed to \code{glm} or \code{glm.nb}
##' @return A \code{SingleSimFit} object
##' @seealso \code{\link{SingleSimFit.object}}
##' @export
Simfit <- function(x,equal.dispersion,formula,...){
  UseMethod("Simfit")
}

##' @export
Simfit.default <- function(x,equal.dispersion,formula,...){
  stop("Invalid x for Simfit")
}

