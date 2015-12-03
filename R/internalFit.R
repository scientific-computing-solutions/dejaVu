# check the arguments to Impute function are valid
validateImputeArguments <- function(fit,impute.mechanism,N){
  
  if(class(fit)!="SingleSimFit"){
    stop("Invalid fit argument, must be of class SingleSimFit")
  }
  
  lapply(fit$impute.parameters,function(x){
    if(is.null(x)) stop("Cannot impute using this SingleSimFit object (a negative binomial model was not fit)")
  })
  
  if(class(impute.mechanism)!="ImputeMechanism"){
    stop("Invalid mpute.mechanism argument, must one of class ImputeMechanism")
  }
  
  if(!is.numeric(N)|| is.na(N) || is.infinite(N) || length(N)>1 || N <= 0 || !.internal.is.wholenumber(N)){
    stop("Invalid argument N, must be positive integer")
  }
  
  
  if(!all(impute.mechanism$cols.needed %in% colnames(fit$singleSim$data))){
    stop("This dropout mechanism requires ",paste(impute.mechanism$cols.needed,collapse=", "),
         "as column names in the simulated data frame")
  }
}