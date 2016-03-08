#helper function for GetgenCoeff
.getPredict <- function(models,data,use.uncertainty,formula){
  data$censored.time <- rep(1,nrow(data))
  
  if(use.uncertainty){
    newcoef1 <- mvrnorm(n=1,mu=coef(models[[1]]),Sigma=vcov(models[[1]]))
    newcoef2 <- if(length(models)==2) mvrnorm(n=1,mu=coef(models[[2]]),Sigma=vcov(models[[2]]))
                else newcoef1
  }
  else{
    newcoef1 <- coef(models[[1]])  
    newcoef2 <- if(length(models)==2) coef(models[[2]])
                else newcoef1
  }
  
  
  if(length(models)==1) models[[2]] <- models[[1]]
  
  data$arm <- factor(rep(0,nrow(data)),levels=c("0","1"))
  d <- model.frame(formula,data)
  
  mu <- exp(model.matrix(formula,d) %*% newcoef1)   
 
  data$arm <- factor(rep(1,nrow(data)),levels=c("0","1"))
  d <- model.frame(formula,data)
  c(mu,exp( model.matrix(formula,d) %*% newcoef2))  
}



#sample gamma uncertainty
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
GetgenCoeff <- function(model,data,equal.dispersion,formula){

  #both arms use same model
  if(equal.dispersion){
    
    #get gamma and sd ready for sampling
    gamma <- model$theta
    gamma.sd <- model$SE.theta/model$theta
    
    
    return(function(use.uncertainty=TRUE){
      gamma <- if(use.uncertainty) rep(getNewGamma(gamma, gamma.sd), 2)
               else rep(gamma,2)   
      mu <- .getPredict(list(model),data,use.uncertainty,formula)
      list(gamma=gamma,
         mu=matrix(mu,ncol = 2,byrow=FALSE))}
    )
  
  }
  
  #Now deal with case each arm is using its own model
  gamma <- vapply(model,function(mod){mod$theta},FUN.VALUE = numeric(1))  
  gamma.sd <- vapply(model,function(mod){mod$SE.theta/mod$theta},FUN.VALUE = numeric(1))
  
  function(use.uncertainty=TRUE){
    if(use.uncertainty){
      gamma <- vapply(seq_along(gamma),function(i) getNewGamma(gamma[i],gamma.sd[i]),FUN.VALUE=numeric(1))
    }   
    mu <- .getPredict(model,data,use.uncertainty,formula)
    list(gamma=gamma,
         mu=matrix(mu,ncol = 2,byrow=FALSE))   
  }
  
  
}