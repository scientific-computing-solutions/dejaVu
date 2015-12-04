##' Scenario object
##' 
##' TODO
##' @name Scenario.object
NULL 


##' @export
as.data.frame.Scenario <- function(x,use.adjusted.pval=FALSE,...){
  
  if(!is.logical(use.adjusted.pval) || length(use.adjusted.pval)>1){
    stop("Invalid argument use.adjusted.pval")
  }
  
  if(use.adjusted.pval && class(x[[1]])=="summary.SingleSimFit") {
    stop("Invalid arguments: cannot use adjusted.pval if multiple imputed test statistic is not used")
  }
  
  .extract <- function(name,fun.val=numeric(1)){
    vapply(x,function(x){x[[name]]},FUN.VALUE = fun.val)
  }
  
  pval.string <-if(use.adjusted.pval)  "adjusted.pval" else "pval"
  
  
  data.frame(replica=1:length(x),
             treatment.effect=.extract("treatment.effect"),
             se=.extract("se"),
             pval=.extract(pval.string),
             theta=.extract("theta")
            )

}



##' Boo
##' @export
CreateScenario <- function(object){

  if(!is.list(object)){
    stop("Invalid argument: object must be a list of either summary.ImputeSimFit or summary.SingleSimFit")
  }    
  
  classes <- lapply(object,class)
  if(!all(classes=="summary.ImputeSimFit") && !(all(classes=="summary.SingleSimFit"))){
    stop("Invalid argument: object must be a list of either summary.ImputeSimFit or summary.SingleSimFit")
  }
  
  class(object) <- "Scenario"                             
  object                             
}