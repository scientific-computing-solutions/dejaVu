#Apply the function f to the list models and repeat the (vector) output
#num.reps times. The ... are additional arguments to f
#return a vector
apply.func <- function(models, f, num.reps, ...){
  rep(unlist(lapply(models,f,...)),num.reps)
}


#helper function for GetgenCoeff
#use.uncertainty is logical argument as to whether the 
#mu values are to include uncertainty or not
#see arguments to GetgenCoeff for further details 
.getPredict <- function(models,data,use.uncertainty,formula){
  
  #set newcoef to be a vector of length 2*#covariates,
  #the first half is for model with arm=0, the second 
  #half is for model with arm=1
  if(use.uncertainty){
    newcoef <- apply.func(models,function(model){
      mvrnorm(n=1,mu=coef(model),Sigma=vcov(model))
    },num.reps=(3-length(models)))
  }
  else{
    newcoef <- apply.func(models,coef,num.reps=(3-length(models)))  
  }
  
  #split newcoef into a list of length 2, each containing coefficients
  #for a single treatment arm
  l <- length(newcoef)
  newcoef <- list(newcoef[1:(l/2)],
                  newcoef[((l/2)+1):l])
  
  
  data$censored.time <- rep(1,nrow(data))
  unlist(lapply(0:1,function(i){
    data$arm <- factor(rep(i,nrow(data)),levels=c("0","1"))
    d <- model.frame(formula,data)
    exp(model.matrix(formula,d) %*% newcoef[[1+i]])   
  }))

}


#sample gamma uncertainty
getNewGamma <- function(gamma, gamma.sd){
  exp(rnorm(1, mean=log(gamma), sd=gamma.sd))
}


#given a list of models  from creating a SimFit object
#output a function which will return a list of impute parameters which contains
#a vector of gamma (1/dispersion) values, control then active arm
#a matrix of mu (mean of negative binomial) values, one row per subject
#first column control arm, second column active arm
#the function has a single argument, use.uncertainty if TRUE then
#the function is stochastic and samples the parameters with uncertainty
GetgenCoeff <- function(model,data,formula){

  #extract the 1/dispersion estimate and its standard error 
  gamma <- apply.func(model, "[[", num.reps = 1, "theta")
  gamma.sd <- apply.func(model, "[[", num.reps = 1, "SE.theta")/gamma
  
  function(use.uncertainty=TRUE){
    
    if(!use.uncertainty){
      warning("Not using uncertainty in parameters therefore imputation will not be proper")
    }
    
    #if we use uncertainty then sample gamma
    if(use.uncertainty){
      gamma <- vapply(seq_along(gamma),
                  function(i) getNewGamma(gamma[i],gamma.sd[i]),FUN.VALUE=numeric(1))
    }
    
    #need to ensure gamma has length 2, one for each arm
    if(length(gamma)==1) gamma <- rep(gamma,2)
    
    #and get the predicted means
    mu <- .getPredict(model,data,use.uncertainty,formula)
    
    list(gamma=gamma, mu=matrix(mu,ncol = 2,byrow=FALSE))   
  }

}