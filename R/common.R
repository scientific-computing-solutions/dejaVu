#cat a named list
.internal.output.list <- function(x){ 
  if(!is.null(x)){
    cat("Parameters:",fill=TRUE)
    invisible(lapply(seq_along(x),function(i){
      cat("  ",names(x)[i],": ",x[[i]],"\n",sep="")}))
  }
}

.internal.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol