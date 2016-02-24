##' @importFrom MASS glm.nb
NULL

#cat a named list
.internal.output.list <- function(x){ 
  if(!is.null(x)){
    cat("Parameters:",fill=TRUE)
    invisible(lapply(seq_along(x),function(i){
      cat("  ",names(x)[i],": ",x[[i]],"\n",sep=" ")}))
  }
}

.internal.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

.internal.is.finite.number <- function(x){
  if(!is.numeric(x) || is.na(x) || length(x)>1 || is.infinite(x)){
    return(FALSE)
  }
  return(TRUE)
}


##' S3 generic to output the number of subjects in each arm for a given object
##' @param x The object
##' @return A vector of the number of subjects
##' @export
subjectsPerArm <- function(x){
  UseMethod("subjectsPerArm")
}

##' @export
subjectsPerArm.default <- function(x){
  stop("Invalid x for calculating number of subjects")
}


##' S3 generic to output the number of subjects in a given object
##' @param x The object
##' @return The number of subjects
##' @export
numberSubjects <- function(x){
  UseMethod("numberSubjects")
}

##' @export
numberSubjects.default <- function(x){
  sum(subjectsPerArm(x))
}

##' S3 generic for fitting models
##' @param x The S3 object
##' @param family Either "negbin" for fitting a negative binomial model (using \code{MASS::glm.nb}),
##' "poisson" for fitting a poisson model (\code{glm}) or "quasipoisson" for fitting a quasipoisson model
##' \code{glm} 
##' @param equal.dispersion logical, should the arms have the same dispersion parameter when
##' fitting negative binomial models
##' @param covar A formula containing the additional covariates to be used when calling \code{glm.nb} if 
##' no covariates are included in the model this should be NULL, for example \code{~covar1 + covar2}
##' See vignette for further details  
##' @param ... Additional arguments to be passed to \code{glm} or \code{glm.nb}
##' @return A \code{SingleSimFit} object
##' @seealso \code{\link{SingleSimFit.object}}
##' @export
Simfit <- function(x,family="negbin",equal.dispersion=TRUE,covar=NULL,...){
  UseMethod("Simfit")
}

##' @export
Simfit.default <- function(x,family="negbin",equal.dispersion=TRUE,covar=NULL,...){
  stop("Invalid x for Simfit")
}


.onAttach <- function(libname,pkgname){
  packageStartupMessage("Warning: This package is in development. DO NOT USE without contacting David Ruau or Paul Metcalfe. ",
                        "Note due to namespace conflicts between MASS and trafficlight this package will not work if trafficlight ",
                        "is loaded")
}
