#helper function for GetgenCoeff
.getPredict <- function(models,data){
  data$censored.time <- rep(1,nrow(data))
  data$arm <- factor(rep(0,nrow(data)))
  mu <- exp(predict(models[[1]],newdata=data,type="link"))
  data$arm <- factor(rep(1,nrow(data)))  
  c(mu,exp(predict(models[[2]],newdata=data,type="link")))  
}


getNewGamma <- function(gamma, gamma.sd){
  exp(rnorm(1, mean=log(gamma), sd=gamma.sd))
}


#given model (or list of models) of from creating a SimFit object
#output a function which will return a list of impute parameters which contains
#a vector of gamma (1/dispersion) values, control then active arm
#a matrix of mu (mean of negative binomial) values, one row per subject
#first column control arm, second column active arm
#the function has a single argument, use.uncertainty if TRUE then
#the function is stochastic and samples the parameters with uncertainty
GetgenCoeff <- function(model,data,equal.dispersion){

  if(equal.dispersion){
    gamma <- model$theta
    gamma.sd <- model$SE.theta/model$theta
    mu<- .getPredict(list(model,model),data)
    return(function(use.uncertainty=TRUE){
      gamma <- if(use.uncertainty) rep(getNewGamma(gamma, gamma.sd), 2)
               else rep(gamma,2)   
      list(gamma=gamma,
         mu=matrix(mu,ncol = 2,byrow=FALSE))}
    )
  }
  
  
  gamma <- vapply(model,function(mod){mod$theta},FUN.VALUE = numeric(1))  
  gamma.sd <- vapply(model,function(mod){mod$SE.theta/mod$theta},FUN.VALUE = numeric(1))
  mu <-  .getPredict(model,data) 
  
  function(use.uncertainty=TRUE){
    if(use.uncertainty){
      gamma <- vapply(seq_along(gamma),function(i) getNewGamma(gamma[i],gamma.sd[i]),FUN.VALUE=numeric(1))
    }   
    
    list(gamma=gamma,
         mu=matrix(mu,ncol = 2,byrow=FALSE))   
  }
  
  
}