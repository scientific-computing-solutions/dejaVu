##' @importFrom MASS glm.nb
NULL

#' @import stats
#' @import utils

##' @importFrom MASS mvrnorm

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
##' @return A vector of the number of subjects in each arm
##' @examples 
##' sim <- SimulateComplete(study.time=365,number.subjects=50,
##' event.rates=c(0.01,0.005),dispersions=0.25)
##' subjectsPerArm(sim)
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
##' @examples 
##' sim <- SimulateComplete(study.time=365,number.subjects=50,
##' event.rates=c(0.01,0.005),dispersions=0.25)
##' numberSubjects(sim)
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
##' 
##' @examples  
##' set.seed(1234)
##' sim <- SimulateComplete(study.time=1,number.subjects=50,
##'       event.rates=c(0.1,0.05),dispersions=0.1)
##' summary(Simfit(sim,equal.dispersion=TRUE))
##'
##' 
Simfit <- function(x,family="negbin",equal.dispersion=TRUE,covar=NULL,...){
  UseMethod("Simfit")
}

##' @export
Simfit.default <- function(x,family="negbin",equal.dispersion=TRUE,covar=NULL,...){
  stop("Invalid x for Simfit")
}


.onAttach <- function(libname,pkgname){
  packageStartupMessage("Note due to namespace conflicts between MASS and trafficlight this package will not work if trafficlight ",
                        "is loaded")
}
