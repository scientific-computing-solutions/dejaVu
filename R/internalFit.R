#If no formula is given to the fit function then 
#a default formula is used (different depending on whether the arms
#have the same dispersion parameters)
GetDefaultFormula <- function(equal.dispersion){
  if(!is.logical(equal.dispersion) || length(equal.dispersion)>1){
    stop("Invalid equal.dispersion argument")
  }
  
  if(equal.dispersion){
    return(formula("observed.events ~ arm + offset(log(censored.time))"))
  }
  
  formula("observed.events ~ offset(log(censored.time))")
  
} 

#' check the arguments to Impute function are valid
validateImputeArguments <- function(fit,impute.mechanism,N){
  
  if(class(fit)!="SingleSimFit"){
    stop("Invalid fit argument, must be of class SingleSimFit")
  }
  
  if(class(impute.mechanism)!="ImputeMechanism"){
    stop("Invalid mpute.mechanism argument, must one of class ImputeMechanism")
  }
  
  if(!is.numeric(N)|| is.na(N) || is.infinite(N) || length(N)>1 || N <= 0 || !.internal.is.wholenumber(N)){
    stop("Invalid argument N, must be positive integer")
  }
  
  
  if(!all(impute.mechanism$cols.needed %in% colnames(fit$singleSim$data))){
    stop("This dropout mechanism requires ",paste(drop.mechanism$cols.needed,collapse=", "),
         "as column names in the simulated data frame")
  }
}