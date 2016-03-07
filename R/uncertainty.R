#helper function for GetGamma_mu
.getPredict <- function(models,data){
  data$censored.time <- rep(1,nrow(data))
  data$arm <- factor(rep(0,nrow(data)))
  mu <- exp(predict(models[[1]],newdata=data,type="link"))
  data$arm <- factor(rep(1,nrow(data)))  
  c(mu,exp(predict(models[[2]],newdata=data,type="link")))  
}

#given model (or list of models) of from creating a SimFit object
#output a list of impute.parameters which contains
#the value of equal.dispersion
#a vector of gamma (1/dispersion) values, control then active arm
#a matrix of mu (mean of negative binomial) values, one row per subject
#first column control arm, second column active arm
GetGamma_mu <- function(model,data,equal.dispersion){
  if(equal.dispersion){
    mod.summary <- summary(model)
    gamma <- rep(mod.summary$theta,2)
    mu<- .getPredict(list(model,model),data)
  }
  else{
    gamma <- vapply(model,function(mod){summary(mod)$theta},FUN.VALUE = numeric(1)) 
    mu <-  .getPredict(model,data) 
  }
  #when this function is changed then don't forget to update SingleSimFit.object in Roxygen
  #and the user guide
  gamma_mu_function <- function(){list(gamma=gamma,mu=matrix(mu,ncol = 2,byrow=FALSE))}
  
  return(list(gamma_mu_function=gamma_mu_function,
              equal.dispersion=equal.dispersion))
}