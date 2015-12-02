GetDefaultFormula <- function(equal.dispersion){
  if(!is.logical(equal.dispersion) || length(equal.dispersion)>1){
    stop("Invalid equal.dispersion argument")
  }
  
  if(equal.dispersion){
    return(formula("observed.events ~ arm + offset(log(censored.time))"))
  }
  
  formula("observed.events ~ offset(log(censored.time))")
  
} 

validateImputeArguments <- function(fit,impute.mechanism,N){
  #TODO
}