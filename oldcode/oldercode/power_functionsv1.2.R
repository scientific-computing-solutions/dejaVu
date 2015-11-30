# Missing Data Project Simulation Program
# 

library(MASS)
library(plyr)

#complete data

complete <- function(n, rate.placebo, rr, T, dispersion.placebo, 
dispersion.treatment, N.rep, alpha){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  


  poisson.se <- rep(NA,N.rep)
  poisson.se2 <- rep(NA,N.rep)
  poisson.rratio <- rep(NA,N.rep)
  poisson.pval <- rep(NA,N.rep)
  quasipoisson.se <- rep(NA,N.rep)
  quasipoisson.se2 <- rep(NA,N.rep)
  quasipoisson.rratio <- rep(NA,N.rep)
  quasipoisson.pval <- rep(NA,N.rep)
  nbd.se <- rep(NA,N.rep)
  nbd.se2 <- rep(NA,N.rep)
  nbd.rratio <- rep(NA,N.rep)
  nbd.pval <- rep(NA,N.rep)
  nbd.theta <- rep(NA,N.rep)
  nbd.warn <- rep(NA,N.rep)

  err <- rep(NA,N.rep)
  
  for (i in 1:N.rep){

    err[i] <-tryCatch({
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- rep(T,n)    
    dropout_time.treatment <- rep(T,n)
    
    dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    poisson.rratio[i] <- exp(poisson.summary$coeffi[2,1])
    poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
    poisson.se2[i] <- poisson.summary$coeffi[2,2]
    poisson.pval[i] <- poisson.summary$coef[2,4] #p-value
    
    quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    quasipoisson.rratio[i] <- exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]
    quasipoisson.pval[i] <- quasipoisson.summary$coef[2,4] 
    
    nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    nbd.warn[i] <- (!is.null(nbd.summary$th.warn))
    nbd.rratio[i] <- exp(nbd.summary$coeffi[2,1])
    nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])
    nbd.se2[i] <- nbd.summary$coeffi[2,2]
    nbd.pval[i] <- nbd.summary$coef[2,4] 
    nbd.theta[i] <- nbd.summary$theta

    err[i] <- 0
    }, warning = function(w){return(1)})
  }


  
  # excludes cases with warnings

  powerpois<-sum(subset(poisson.pval,err==0)<alpha)/length(subset(poisson.pval,err==0))
  powerqpois<-sum(subset(quasipoisson.pval,err==0)<alpha)/length(subset(quasipoisson.pval,err==0))
  powernb<-sum(subset(nbd.pval,err==0)<alpha)/length(subset(nbd.pval,err==0))
  power <- as.data.frame(cbind(powerpois,powerqpois,powernb))
  
  rratio <- as.data.frame(cbind(poisson.rratio,quasipoisson.rratio,nbd.rratio))
  se <- as.data.frame(cbind(poisson.se,quasipoisson.se,nbd.se,poisson.se2,quasipoisson.se2,nbd.se2))
  
  quannbtheta<-quantile(nbd.theta)
  nbthetawarn <- sum(nbd.warn)
  
  quandrop.p <- quantile(dropout.p)
  quandrop.t <- quantile(dropout.t)
  dropout <- as.data.frame(cbind(quandrop.p,quandrop.t))
  
#  result <- list(err=err,power=power,theta=quannbtheta,rratio=rratio,se=se,nbthetawarn=nbthetawarn,dropout=dropout)
     result <- list(err=err,power=power,poisson.pval=poisson.pval,quasipoisson.pval=quasipoisson.pval,nbd.pval=nbd.pval,theta=quannbtheta,rratio=rratio,se=se,nbthetawarn=nbthetawarn,dropout=dropout)
  return(result)
  
}

#
# MCAR without variation in drop-out rate
# 
#
# 
sim_mcar1<-function(n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,drop_rate.placebo,drop_rate.treatment,N.rep,N.mi,alpha){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  censor.time.p <- rep(NA,N.rep)
  censor.time.t <- rep(NA,N.rep)
  ERgt75p.p <- rep(NA,N.rep)
  ERgt75p.t <- rep(NA,N.rep)
  ERle75p.p <- rep(NA,N.rep)
  ERle75p.t <- rep(NA,N.rep)
  
  # no imputation (direct likelihood approach)  
  dl_poisson.se <- rep(NA,N.rep)
  dl_poisson.rratio <- rep(NA,N.rep)
  dl_poisson.pval <- rep(NA,N.rep)
  
  dl_quasipoisson.se <- rep(NA,N.rep)
  dl_quasipoisson.rratio <- rep(NA,N.rep)
  dl_quasipoisson.pval <- rep(NA,N.rep)
  
  
  dl_nbd.se <- rep(NA,N.rep)
  dl_nbd.rratio <- rep(NA,N.rep)
  dl_nbd.pval <- rep(NA,N.rep)
  dl_nbd.theta <- rep(NA,N.rep)
  dl_nbd.warn <- rep(NA,N.rep)
  
  # jump to reference 
  j2r_poisson.se <- rep(NA,N.rep)
  j2r_poisson.rratio <- rep(NA,N.rep)
  j2r_poisson.pval <- rep(NA,N.rep)
  j2r_poisson.apval <- rep(NA,N.rep)
  j2r_poisson.df <- rep(NA,N.rep)
  j2r_poisson.adf <- rep(NA,N.rep) 
  
  j2r_quasipoisson.se <- rep(NA,N.rep)
  j2r_quasipoisson.rratio <- rep(NA,N.rep)
  j2r_quasipoisson.pval <- rep(NA,N.rep)
  j2r_quasipoisson.apval <- rep(NA,N.rep)
  j2r_quasipoisson.df <- rep(NA,N.rep)
  j2r_quasipoisson.adf <- rep(NA,N.rep)
  
  j2r_nbd.se <- rep(NA,N.rep)
  j2r_nbd.rratio <- rep(NA,N.rep)
  j2r_nbd.pval <- rep(NA,N.rep)
  j2r_nbd.apval <- rep(NA,N.rep)
  j2r_nbd.df <- rep(NA,N.rep)
  j2r_nbd.adf <- rep(NA,N.rep)
  
  # efficacy estimand
  efy_poisson.se <- rep(NA,N.rep)
  efy_poisson.rratio <- rep(NA,N.rep)
  efy_poisson.pval <- rep(NA,N.rep)
  efy_poisson.apval <- rep(NA,N.rep)
  efy_poisson.df <- rep(NA,N.rep)
  efy_poisson.adf <- rep(NA,N.rep)
  
  efy_quasipoisson.se <- rep(NA,N.rep)
  efy_quasipoisson.rratio <- rep(NA,N.rep)
  efy_quasipoisson.pval <- rep(NA,N.rep)
  efy_quasipoisson.apval <- rep(NA,N.rep)
  efy_quasipoisson.df <- rep(NA,N.rep)
  efy_quasipoisson.adf <- rep(NA,N.rep)
  
  efy_nbd.se <- rep(NA,N.rep)
  efy_nbd.rratio <- rep(NA,N.rep)
  efy_nbd.pval <- rep(NA,N.rep)
  efy_nbd.apval <- rep(NA,N.rep)
  efy_nbd.df <- rep(NA,N.rep)
  efy_nbd.adf <- rep(NA,N.rep)
  
  
  for (i in 1:N.rep){
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- dropout_time1(T,n,drop_rate.placebo)
    dropout_time.treatment <- dropout_time1(T,n,drop_rate.treatment)
    
    dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    # mean ER of patients who stayed for at least 75% of total follow-up period vs early dropouts
    ERgt75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]>0.75*T))
    ERgt75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]>0.75*T))
    ERle75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<=0.75*T))
    ERle75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<=0.75*T))
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p[i] <- mean(subset(nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<T))
    censor.time.t[i] <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    # no imputation (direct likelihood approach)
    
    dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]
    dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
    dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
    
    dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]
    dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
    dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
    
    dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
    dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]
    dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
    dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
    dl_nbd.theta[i] <- dl_nbd.summary$theta
    
    # jump to reference
    
    j2r <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    j2r_poisson.se[i] <- j2r$mi_poisson.se
    j2r_poisson.rratio[i] <- j2r$mi_poisson.trt
    j2r_poisson.df[i] <- j2r$mi_poisson.df
    j2r_poisson.adf[i] <- j2r$mi_poisson.adf
    j2r_poisson.pval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se[i]),df=j2r_poisson.df[i]))
    j2r_poisson.apval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se[i]),df=j2r_poisson.adf[i]))
    
    j2r_quasipoisson.se[i] <- j2r$mi_quasipoisson.se
    j2r_quasipoisson.rratio[i] <- j2r$mi_quasipoisson.trt
    j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df
    j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf
    j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se[i]),df=j2r_quasipoisson.df[i]))
    j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se[i]),df=j2r_quasipoisson.adf[i]))
    
    j2r_nbd.se[i] <- j2r$mi_nbd.se
    j2r_nbd.rratio[i] <- j2r$mi_nbd.trt
    j2r_nbd.df[i] <- j2r$mi_nbd.df
    j2r_nbd.adf[i] <- j2r$mi_nbd.adf
    j2r_nbd.pval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se[i]),df=j2r_nbd.df[i]))
    j2r_nbd.apval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se[i]),df=j2r_nbd.adf[i]))
    
    # efficacy
    
    efy <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)
    efy_poisson.se[i] <- efy$mi_poisson.se
    efy_poisson.rratio[i] <- efy$mi_poisson.trt
    efy_poisson.df[i] <- efy$mi_poisson.df
    efy_poisson.adf[i] <- efy$mi_poisson.adf
    efy_poisson.pval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se[i]),df=efy_poisson.df[i]))
    efy_poisson.apval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se[i]),df=efy_poisson.adf[i]))
    
    efy_quasipoisson.se[i] <- efy$mi_quasipoisson.se
    efy_quasipoisson.rratio[i] <- efy$mi_quasipoisson.trt
    efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df
    efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf
    efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se[i]),df=efy_quasipoisson.df[i]))
    efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se[i]),df=efy_quasipoisson.adf[i]))
    
    efy_nbd.se[i] <- efy$mi_nbd.se
    efy_nbd.rratio[i] <- efy$mi_nbd.trt
    efy_nbd.df[i] <- efy$mi_nbd.df
    efy_nbd.adf[i] <- efy$mi_nbd.adf
    efy_nbd.pval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se[i]),df=efy_nbd.df[i]))
    efy_nbd.apval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se[i]),df=efy_nbd.adf[i]))
    
    
    
  }
  
  dl_powerpois<-sum(dl_poisson.pval<alpha)/N.rep 
  dl_powerqpois<-sum(dl_quasipoisson.pval<alpha)/N.rep 
  dl_powernb<-sum(dl_nbd.pval<alpha)/N.rep 
  
  j2r_powerpois<-sum(j2r_poisson.pval<alpha)/N.rep 
  j2r_powerqpois<-sum(j2r_quasipoisson.pval<alpha)/N.rep 
  j2r_powernb<-sum(j2r_nbd.pval<alpha)/N.rep 
  
  j2r_powerpoisa<-sum(j2r_poisson.apval<alpha)/N.rep 
  j2r_powerqpoisa<-sum(j2r_quasipoisson.apval<alpha)/N.rep 
  j2r_powernba<-sum(j2r_nbd.apval<alpha)/N.rep 
  
  efy_powerpois<-sum(efy_poisson.pval<alpha)/N.rep 
  efy_powerqpois<-sum(efy_quasipoisson.pval<alpha)/N.rep 
  efy_powernb<-sum(efy_nbd.pval<alpha)/N.rep 
  
  efy_powerpoisa<-sum(efy_poisson.apval<alpha)/N.rep 
  efy_powerqpoisa<-sum(efy_quasipoisson.apval<alpha)/N.rep 
  efy_powernba<-sum(efy_nbd.apval<alpha)/N.rep 
  
  power <- as.data.frame(cbind(dl_powerpois,dl_powerqpois,dl_powernb,j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba))
  
  rratio <- as.data.frame(cbind(dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio))
  se <- as.data.frame(cbind(dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,efy_poisson.se,efy_quasipoisson.se,efy_nbd.se))
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf))
  
  dl_summnbtheta<-summary(dl_nbd.theta)
  dl_nbthetawarn <- sum(dl_nbd.warn)/N.rep
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
  
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  result <- list(power=power,dl_theta=dl_summnbtheta,rratio=rratio,se=se,mi_df=mi_df,dl_nbthetawarn=dl_nbthetawarn,dropout=dropout,censortime=censortime,ERgt75p=ERgt75p,ERle75p=ERle75p)
  return(result)
  
  
}

# MCAR with variation in drop-out rates

sim_mcar2<-function(n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,drop_rate.placebo,drop_rate.treatment,drop_rate_variance,N.rep,N.mi,alpha){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  censor.time.p <- rep(NA,N.rep)
  censor.time.t <- rep(NA,N.rep)
  ERgt75p.p <- rep(NA,N.rep)
  ERgt75p.t <- rep(NA,N.rep)
  ERle75p.p <- rep(NA,N.rep)
  ERle75p.t <- rep(NA,N.rep)
  
  # no imputation (direct likelihood approach)
  dl_poisson.se <- rep(NA,N.rep)
  dl_poisson.rratio <- rep(NA,N.rep)
  dl_poisson.pval <- rep(NA,N.rep)
  
  dl_quasipoisson.se <- rep(NA,N.rep)
  dl_quasipoisson.rratio <- rep(NA,N.rep)
  dl_quasipoisson.pval <- rep(NA,N.rep)
  
  dl_nbd.se <- rep(NA,N.rep)
  dl_nbd.rratio <- rep(NA,N.rep)
  dl_nbd.pval <- rep(NA,N.rep)
  dl_nbd.theta <- rep(NA,N.rep)
  dl_nbd.warn <- rep(NA,N.rep)
  
  # jump to reference 
  j2r_poisson.se <- rep(NA,N.rep)
  j2r_poisson.rratio <- rep(NA,N.rep)
  j2r_poisson.pval <- rep(NA,N.rep)
  j2r_poisson.apval <- rep(NA,N.rep)
  j2r_poisson.df <- rep(NA,N.rep)
  j2r_poisson.adf <- rep(NA,N.rep) 
  
  j2r_quasipoisson.se <- rep(NA,N.rep)
  j2r_quasipoisson.rratio <- rep(NA,N.rep)
  j2r_quasipoisson.pval <- rep(NA,N.rep)
  j2r_quasipoisson.apval <- rep(NA,N.rep)
  j2r_quasipoisson.df <- rep(NA,N.rep)
  j2r_quasipoisson.adf <- rep(NA,N.rep)
  
  j2r_nbd.se <- rep(NA,N.rep)
  j2r_nbd.rratio <- rep(NA,N.rep)
  j2r_nbd.pval <- rep(NA,N.rep)
  j2r_nbd.apval <- rep(NA,N.rep)
  j2r_nbd.df <- rep(NA,N.rep)
  j2r_nbd.adf <- rep(NA,N.rep)
  
  # efficacy estimand
  efy_poisson.se <- rep(NA,N.rep)
  efy_poisson.rratio <- rep(NA,N.rep)
  efy_poisson.pval <- rep(NA,N.rep)
  efy_poisson.apval <- rep(NA,N.rep)
  efy_poisson.df <- rep(NA,N.rep)
  efy_poisson.adf <- rep(NA,N.rep)
  
  efy_quasipoisson.se <- rep(NA,N.rep)
  efy_quasipoisson.rratio <- rep(NA,N.rep)
  efy_quasipoisson.pval <- rep(NA,N.rep)
  efy_quasipoisson.apval <- rep(NA,N.rep)
  efy_quasipoisson.df <- rep(NA,N.rep)
  efy_quasipoisson.adf <- rep(NA,N.rep)
  
  efy_nbd.se <- rep(NA,N.rep)
  efy_nbd.rratio <- rep(NA,N.rep)
  efy_nbd.pval <- rep(NA,N.rep)
  efy_nbd.apval <- rep(NA,N.rep)
  efy_nbd.df <- rep(NA,N.rep)
  efy_nbd.adf <- rep(NA,N.rep)
  
  
  for (i in 1:N.rep){
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- dropout_time2(T,n,drop_rate.placebo,drop_rate_variance)
    dropout_time.treatment <- dropout_time2(T,n,drop_rate.treatment,drop_rate_variance)
    
    dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    # mean ER of patients who stayed for at least 75% of total follow-up period vs early dropouts
    ERgt75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]>0.75*T))
    ERgt75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]>0.75*T))
    ERle75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<=0.75*T))
    ERle75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<=0.75*T))
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p[i] <- mean(subset(nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<T))
    censor.time.t[i] <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    # no imputation (direct likelihood approach)
    
    dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]
    dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
    dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
    
    dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]
    dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
    dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
    
    dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
    dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]
    dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
    dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
    dl_nbd.theta[i] <- dl_nbd.summary$theta
    
    # jump to reference
    
    j2r <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    j2r_poisson.se[i] <- j2r$mi_poisson.se
    j2r_poisson.rratio[i] <- j2r$mi_poisson.trt
    j2r_poisson.df[i] <- j2r$mi_poisson.df
    j2r_poisson.adf[i] <- j2r$mi_poisson.adf
    j2r_poisson.pval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se[i]),df=j2r_poisson.df[i]))
    j2r_poisson.apval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se[i]),df=j2r_poisson.adf[i]))
    
    j2r_quasipoisson.se[i] <- j2r$mi_quasipoisson.se
    j2r_quasipoisson.rratio[i] <- j2r$mi_quasipoisson.trt
    j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df
    j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf
    j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se[i]),df=j2r_quasipoisson.df[i]))
    j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se[i]),df=j2r_quasipoisson.adf[i]))
    
    j2r_nbd.se[i] <- j2r$mi_nbd.se
    j2r_nbd.rratio[i] <- j2r$mi_nbd.trt
    j2r_nbd.df[i] <- j2r$mi_nbd.df
    j2r_nbd.adf[i] <- j2r$mi_nbd.adf
    j2r_nbd.pval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se[i]),df=j2r_nbd.df[i]))
    j2r_nbd.apval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se[i]),df=j2r_nbd.adf[i]))
    
    # efficacy
    
    efy <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)
    efy_poisson.se[i] <- efy$mi_poisson.se
    efy_poisson.rratio[i] <- efy$mi_poisson.trt
    efy_poisson.df[i] <- efy$mi_poisson.df
    efy_poisson.adf[i] <- efy$mi_poisson.adf
    efy_poisson.pval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se[i]),df=efy_poisson.df[i]))
    efy_poisson.apval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se[i]),df=efy_poisson.adf[i]))
    
    efy_quasipoisson.se[i] <- efy$mi_quasipoisson.se
    efy_quasipoisson.rratio[i] <- efy$mi_quasipoisson.trt
    efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df
    efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf
    efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se[i]),df=efy_quasipoisson.df[i]))
    efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se[i]),df=efy_quasipoisson.adf[i]))
    
    efy_nbd.se[i] <- efy$mi_nbd.se
    efy_nbd.rratio[i] <- efy$mi_nbd.trt
    efy_nbd.df[i] <- efy$mi_nbd.df
    efy_nbd.adf[i] <- efy$mi_nbd.adf
    efy_nbd.pval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se[i]),df=efy_nbd.df[i]))
    efy_nbd.apval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se[i]),df=efy_nbd.adf[i]))
    
    
    
  }
  
  dl_powerpois<-sum(dl_poisson.pval<alpha)/N.rep 
  dl_powerqpois<-sum(dl_quasipoisson.pval<alpha)/N.rep 
  dl_powernb<-sum(dl_nbd.pval<alpha)/N.rep 
  
  j2r_powerpois<-sum(j2r_poisson.pval<alpha)/N.rep 
  j2r_powerqpois<-sum(j2r_quasipoisson.pval<alpha)/N.rep 
  j2r_powernb<-sum(j2r_nbd.pval<alpha)/N.rep 
  
  j2r_powerpoisa<-sum(j2r_poisson.apval<alpha)/N.rep 
  j2r_powerqpoisa<-sum(j2r_quasipoisson.apval<alpha)/N.rep 
  j2r_powernba<-sum(j2r_nbd.apval<alpha)/N.rep 
  
  efy_powerpois<-sum(efy_poisson.pval<alpha)/N.rep 
  efy_powerqpois<-sum(efy_quasipoisson.pval<alpha)/N.rep 
  efy_powernb<-sum(efy_nbd.pval<alpha)/N.rep 
  
  efy_powerpoisa<-sum(efy_poisson.apval<alpha)/N.rep 
  efy_powerqpoisa<-sum(efy_quasipoisson.apval<alpha)/N.rep 
  efy_powernba<-sum(efy_nbd.apval<alpha)/N.rep 
  
  power <- as.data.frame(cbind(dl_powerpois,dl_powerqpois,dl_powernb,j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba))
  
  rratio <- as.data.frame(cbind(dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio))
  se <- as.data.frame(cbind(dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,efy_poisson.se,efy_quasipoisson.se,efy_nbd.se))
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf))
  
  dl_summnbtheta<-summary(dl_nbd.theta)
  dl_nbthetawarn <- sum(dl_nbd.warn)/N.rep
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
  
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  result <- list(power=power,dl_theta=dl_summnbtheta,rratio=rratio,se=se,mi_df=mi_df,dl_nbthetawarn=dl_nbthetawarn,dropout=dropout,censortime=censortime,ERgt75p=ERgt75p,ERle75p=ERle75p)
  return(result)
  
  
}


# MAR without variation in drop-out rate
# 

sim_mar1<-function(n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,N.rep,N.mi,alpha){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  censor.time.p <- rep(NA,N.rep)
  censor.time.t <- rep(NA,N.rep)
  ERgt75p.p <- rep(NA,N.rep)
  ERgt75p.t <- rep(NA,N.rep)
  ERle75p.p <- rep(NA,N.rep)
  ERle75p.t <- rep(NA,N.rep)
  
  
  # no imputation (direct likelihood approach)
  dl_poisson.se <- rep(NA,N.rep)
  dl_poisson.rratio <- rep(NA,N.rep)
  dl_poisson.pval <- rep(NA,N.rep)
  
  dl_quasipoisson.se <- rep(NA,N.rep)
  dl_quasipoisson.rratio <- rep(NA,N.rep)
  dl_quasipoisson.pval <- rep(NA,N.rep)
  
  dl_nbd.se <- rep(NA,N.rep)
  dl_nbd.rratio <- rep(NA,N.rep)
  dl_nbd.pval <- rep(NA,N.rep)
  dl_nbd.theta <- rep(NA,N.rep)
  dl_nbd.warn <- rep(NA,N.rep)
  
  # jump to reference 
  j2r_poisson.se <- rep(NA,N.rep)
  j2r_poisson.rratio <- rep(NA,N.rep)
  j2r_poisson.pval <- rep(NA,N.rep)
  j2r_poisson.apval <- rep(NA,N.rep)
  j2r_poisson.df <- rep(NA,N.rep)
  j2r_poisson.adf <- rep(NA,N.rep) 
  
  j2r_quasipoisson.se <- rep(NA,N.rep)
  j2r_quasipoisson.rratio <- rep(NA,N.rep)
  j2r_quasipoisson.pval <- rep(NA,N.rep)
  j2r_quasipoisson.apval <- rep(NA,N.rep)
  j2r_quasipoisson.df <- rep(NA,N.rep)
  j2r_quasipoisson.adf <- rep(NA,N.rep)
  
  j2r_nbd.se <- rep(NA,N.rep)
  j2r_nbd.rratio <- rep(NA,N.rep)
  j2r_nbd.pval <- rep(NA,N.rep)
  j2r_nbd.apval <- rep(NA,N.rep)
  j2r_nbd.df <- rep(NA,N.rep)
  j2r_nbd.adf <- rep(NA,N.rep)
  
  # efficacy estimand
  efy_poisson.se <- rep(NA,N.rep)
  efy_poisson.rratio <- rep(NA,N.rep)
  efy_poisson.pval <- rep(NA,N.rep)
  efy_poisson.apval <- rep(NA,N.rep)
  efy_poisson.df <- rep(NA,N.rep)
  efy_poisson.adf <- rep(NA,N.rep)
  
  efy_quasipoisson.se <- rep(NA,N.rep)
  efy_quasipoisson.rratio <- rep(NA,N.rep)
  efy_quasipoisson.pval <- rep(NA,N.rep)
  efy_quasipoisson.apval <- rep(NA,N.rep)
  efy_quasipoisson.df <- rep(NA,N.rep)
  efy_quasipoisson.adf <- rep(NA,N.rep)
  
  efy_nbd.se <- rep(NA,N.rep)
  efy_nbd.rratio <- rep(NA,N.rep)
  efy_nbd.pval <- rep(NA,N.rep)
  efy_nbd.apval <- rep(NA,N.rep)
  efy_nbd.df <- rep(NA,N.rep)
  efy_nbd.adf <- rep(NA,N.rep)
  
  
  for (i in 1:N.rep){
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- dropout_time3(T,n,nb_complete.placebo,rates)
    dropout_time.treatment <- dropout_time3(T,n,nb_complete.treatment,rates)
    
    dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    # mean ER of patients who stayed for at least 75% of total follow-up period vs early dropouts
    ERgt75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]>0.75*T))
    ERgt75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]>0.75*T))
    ERle75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<=0.75*T))
    ERle75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<=0.75*T))
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p[i] <- mean(subset(nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<T))
    censor.time.t[i] <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    # no imputation (direct likelihood approach)
    
    dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]
    dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
    dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
    
    dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]
    dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
    dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
    
    dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
    dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]
    dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
    dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
    dl_nbd.theta[i] <- dl_nbd.summary$theta
    
    # jump to reference
    
    j2r <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    j2r_poisson.se[i] <- j2r$mi_poisson.se
    j2r_poisson.rratio[i] <- j2r$mi_poisson.trt
    j2r_poisson.df[i] <- j2r$mi_poisson.df
    j2r_poisson.adf[i] <- j2r$mi_poisson.adf
    j2r_poisson.pval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se[i]),df=j2r_poisson.df[i]))
    j2r_poisson.apval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se[i]),df=j2r_poisson.adf[i]))
    
    j2r_quasipoisson.se[i] <- j2r$mi_quasipoisson.se
    j2r_quasipoisson.rratio[i] <- j2r$mi_quasipoisson.trt
    j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df
    j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf
    j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se[i]),df=j2r_quasipoisson.df[i]))
    j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se[i]),df=j2r_quasipoisson.adf[i]))
    
    j2r_nbd.se[i] <- j2r$mi_nbd.se
    j2r_nbd.rratio[i] <- j2r$mi_nbd.trt
    j2r_nbd.df[i] <- j2r$mi_nbd.df
    j2r_nbd.adf[i] <- j2r$mi_nbd.adf
    j2r_nbd.pval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se[i]),df=j2r_nbd.df[i]))
    j2r_nbd.apval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se[i]),df=j2r_nbd.adf[i]))
    
    # efficacy
    
    efy <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)
    efy_poisson.se[i] <- efy$mi_poisson.se
    efy_poisson.rratio[i] <- efy$mi_poisson.trt
    efy_poisson.df[i] <- efy$mi_poisson.df
    efy_poisson.adf[i] <- efy$mi_poisson.adf
    efy_poisson.pval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se[i]),df=efy_poisson.df[i]))
    efy_poisson.apval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se[i]),df=efy_poisson.adf[i]))
    
    efy_quasipoisson.se[i] <- efy$mi_quasipoisson.se
    efy_quasipoisson.rratio[i] <- efy$mi_quasipoisson.trt
    efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df
    efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf
    efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se[i]),df=efy_quasipoisson.df[i]))
    efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se[i]),df=efy_quasipoisson.adf[i]))
    
    efy_nbd.se[i] <- efy$mi_nbd.se
    efy_nbd.rratio[i] <- efy$mi_nbd.trt
    efy_nbd.df[i] <- efy$mi_nbd.df
    efy_nbd.adf[i] <- efy$mi_nbd.adf
    efy_nbd.pval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se[i]),df=efy_nbd.df[i]))
    efy_nbd.apval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se[i]),df=efy_nbd.adf[i]))
    
    
    
  }
  
  dl_powerpois<-sum(dl_poisson.pval<alpha)/N.rep 
  dl_powerqpois<-sum(dl_quasipoisson.pval<alpha)/N.rep 
  dl_powernb<-sum(dl_nbd.pval<alpha)/N.rep 
  
  j2r_powerpois<-sum(j2r_poisson.pval<alpha)/N.rep 
  j2r_powerqpois<-sum(j2r_quasipoisson.pval<alpha)/N.rep 
  j2r_powernb<-sum(j2r_nbd.pval<alpha)/N.rep 
  
  j2r_powerpoisa<-sum(j2r_poisson.apval<alpha)/N.rep 
  j2r_powerqpoisa<-sum(j2r_quasipoisson.apval<alpha)/N.rep 
  j2r_powernba<-sum(j2r_nbd.apval<alpha)/N.rep 
  
  efy_powerpois<-sum(efy_poisson.pval<alpha)/N.rep 
  efy_powerqpois<-sum(efy_quasipoisson.pval<alpha)/N.rep 
  efy_powernb<-sum(efy_nbd.pval<alpha)/N.rep 
  
  efy_powerpoisa<-sum(efy_poisson.apval<alpha)/N.rep 
  efy_powerqpoisa<-sum(efy_quasipoisson.apval<alpha)/N.rep 
  efy_powernba<-sum(efy_nbd.apval<alpha)/N.rep 
  
  power <- as.data.frame(cbind(dl_powerpois,dl_powerqpois,dl_powernb,j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba))
  
  rratio <- as.data.frame(cbind(dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio))
  se <- as.data.frame(cbind(dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,efy_poisson.se,efy_quasipoisson.se,efy_nbd.se))
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf))
  
  dl_summnbtheta<-summary(dl_nbd.theta)
  dl_nbthetawarn <- sum(dl_nbd.warn)/N.rep
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
  
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  result <- list(power=power,dl_theta=dl_summnbtheta,rratio=rratio,se=se,mi_df=mi_df,dl_nbthetawarn=dl_nbthetawarn,dropout=dropout,censortime=censortime,ERgt75p=ERgt75p,ERle75p=ERle75p)
  return(result)
  
  
}






# MAR with variation in drop-out rate
# 


sim_mar2<-function(n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,drop_rate_variance,N.rep,N.mi, alpha, datname){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  
  py.p <- rep(NA,N.rep)
  py.t <- rep(NA,N.rep)
  
  events.dl.p <- rep(NA,N.rep)
  mr.dl.p <- rep(NA,N.rep)
  mr.ef.p <- rep(NA,N.rep)
  events.ef.p <- rep(NA,N.rep)
  mr.co.p <- rep(NA,N.rep)
  events.co.p <- rep(NA,N.rep)
  mr.j2r.p <- rep(NA,N.rep)
  events.j2r.p <- rep(NA,N.rep)
  
  
  events.dl.t <- rep(NA,N.rep)
  mr.dl.t <- rep(NA,N.rep)
  mr.ef.t <- rep(NA,N.rep)
  events.ef.t <- rep(NA,N.rep)
  mr.co.t <- rep(NA,N.rep)
  events.co.t <- rep(NA,N.rep)
  mr.j2r.t <- rep(NA,N.rep)
  events.j2r.t <- rep(NA,N.rep)
  
  censor.time.p <- rep(NA,N.rep)
  censor.time.t <- rep(NA,N.rep)
  t.rate <- rep(NA,N.rep)
  p.rate <- rep(NA,N.rep)
  gamma.t <- rep(NA,N.rep)
  gamma.p <- rep(NA,N.rep)
  
  t.rate2 <- rep(NA,N.rep)
  p.rate2 <- rep(NA,N.rep)
  gamma <- rep(NA,N.rep)
  
  ERgt75p.p <- rep(NA,N.rep)
  ERgt75p.t <- rep(NA,N.rep)
  ERle75p.p <- rep(NA,N.rep)
  ERle75p.t <- rep(NA,N.rep)
  
  #Complete data
  co_poisson.se <- rep(NA,N.rep)
  co_poisson.se2 <- rep(NA,N.rep)
  co_poisson.rratio <- rep(NA,N.rep)
  co_poisson.pval <- rep(NA,N.rep)
  
  co_quasipoisson.se <- rep(NA,N.rep)
  co_quasipoisson.se2 <- rep(NA,N.rep)
  co_quasipoisson.rratio <- rep(NA,N.rep)
  co_quasipoisson.pval <- rep(NA,N.rep)
  
  co_nbd.se <- rep(NA,N.rep)
  co_nbd.se2 <- rep(NA,N.rep)
  co_nbd.rratio <- rep(NA,N.rep)
  co_nbd.pval <- rep(NA,N.rep)
  co_nbd.theta <- rep(NA,N.rep)
  co_nbd.warn <- rep(NA,N.rep)
  
  # no imputation (direct likelihood approach)
  dl_poisson.se <- rep(NA,N.rep)
  dl_poisson.se2 <- rep(NA,N.rep)
  dl_poisson.rratio <- rep(NA,N.rep)
  dl_poisson.pval <- rep(NA,N.rep)
  
  dl_quasipoisson.se <- rep(NA,N.rep)
  dl_quasipoisson.se2 <- rep(NA,N.rep)
  dl_quasipoisson.rratio <- rep(NA,N.rep)
  dl_quasipoisson.pval <- rep(NA,N.rep)
  
  dl_nbd.se <- rep(NA,N.rep)
  dl_nbd.se2 <- rep(NA,N.rep)
  dl_nbd.rratio <- rep(NA,N.rep)
  dl_nbd.pval <- rep(NA,N.rep)
  dl_nbd.theta <- rep(NA,N.rep)
  dl_nbd.warn <- rep(NA,N.rep)
  nbd.warn.p <- rep(NA,N.rep)
  nbd.warn.t <- rep(NA,N.rep)
  
  # jump to reference 
  j2r_poisson.se <- rep(NA,N.rep)
  j2r_poisson.se2 <- rep(NA,N.rep)
  j2r_poisson.rratio <- rep(NA,N.rep)
  j2r_poisson.pval <- rep(NA,N.rep)
  j2r_poisson.apval <- rep(NA,N.rep)
  j2r_poisson.df <- rep(NA,N.rep)
  j2r_poisson.adf <- rep(NA,N.rep) 
  
  j2r_quasipoisson.se <- rep(NA,N.rep)
  j2r_quasipoisson.se2 <- rep(NA,N.rep)
  j2r_quasipoisson.rratio <- rep(NA,N.rep)
  j2r_quasipoisson.pval <- rep(NA,N.rep)
  j2r_quasipoisson.apval <- rep(NA,N.rep)
  j2r_quasipoisson.df <- rep(NA,N.rep)
  j2r_quasipoisson.adf <- rep(NA,N.rep)
  
  j2r_nbd.se <- rep(NA,N.rep)
  j2r_nbd.se2 <- rep(NA,N.rep)
  j2r_nbd.rratio <- rep(NA,N.rep)
  j2r_nbd.pval <- rep(NA,N.rep)
  j2r_nbd.apval <- rep(NA,N.rep)
  j2r_nbd.df <- rep(NA,N.rep)
  j2r_nbd.adf <- rep(NA,N.rep)
  outliers.j2r=c()
  
  # efficacy estimand
  efy_poisson.se <- rep(NA,N.rep)
  efy_poisson.se2 <- rep(NA,N.rep)
  efy_poisson.rratio <- rep(NA,N.rep)
  efy_poisson.pval <- rep(NA,N.rep)
  efy_poisson.apval <- rep(NA,N.rep)
  efy_poisson.df <- rep(NA,N.rep)
  efy_poisson.adf <- rep(NA,N.rep)
  
  efy_quasipoisson.se <- rep(NA,N.rep)
  efy_quasipoisson.se2 <- rep(NA,N.rep)
  efy_quasipoisson.rratio <- rep(NA,N.rep)
  efy_quasipoisson.pval <- rep(NA,N.rep)
  efy_quasipoisson.apval <- rep(NA,N.rep)
  efy_quasipoisson.df <- rep(NA,N.rep)
  efy_quasipoisson.adf <- rep(NA,N.rep)
  
  efy_nbd.se <- rep(NA,N.rep)
  efy_nbd.se2 <- rep(NA,N.rep)
  efy_nbd.rratio <- rep(NA,N.rep)
  efy_nbd.pval <- rep(NA,N.rep)
  efy_nbd.apval <- rep(NA,N.rep)
  efy_nbd.df <- rep(NA,N.rep)
  efy_nbd.adf <- rep(NA,N.rep)
  outliers.ef=c()
  alldata=data.frame()
  
  err <- rep(NA,N.rep)
  
  for (i in 1:N.rep){
    if ((i/10 - round(i/10)) == 0) print(i)
    err[i] <-tryCatch({
    nb_complete.placebo   <- nb_complete(T, n, rate.placebo,   dispersion.placebo)
    nb_complete.treatment <- nb_complete(T, n, rate.treatment, dispersion.treatment)
    
    dropout_time.placebo <- dropout_time4(T,n,nb_complete.placebo,rates,drop_rate_variance)
    dropout_time.treatment <- dropout_time4(T,n,nb_complete.treatment,rates,drop_rate_variance)
    
    dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    cc.p=rep(0,n)
    for (k in 1:n) cc.p[k]=sum(!is.na(nb_complete.placebo[[k]]))
    cc.t=rep(0,n)
    for (k in 1:n) cc.t[k]=sum(!is.na(nb_complete.treatment[[k]]))
    
        
    # mean ER of patients who stayed for at least 75% of total follow-up period vs early dropouts
    ERgt75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]>0.75*T))
    ERgt75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]>0.75*T))
    ERle75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<=0.75*T))
    ERle75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<=0.75*T))
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p[i] <- mean(subset(nb_incomplete.placebo[[2]],  nb_incomplete.placebo[[2]]<T))
    censor.time.t[i] <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    id=1:n
    Xchar=rep("",length(X)) 
    Xchar[X==0]="Placebo"
    Xchar[X==1]="Active"
        
    
    
    #MG add
    py.p[i]=sum(nb_incomplete.placebo[[2]])
    py.t[i]=sum(nb_incomplete.treatment[[2]])
    
    events.co.p[i]=0
    for (k in 1:n) events.co.p[i]=events.co.p[i] + sum(!is.na(nb_complete.placebo[[k]]))                 
    mr.co.p[i]=events.co.p[i]/(n)
    
    events.co.t[i]=0
    for (k in 1:n) events.co.t[i]=events.co.t[i] + sum(!is.na(nb_complete.treatment[[k]]))                
    mr.co.t[i]=events.co.t[i]/(n)
    
    #for (k in 1:n) events.dl.p[i]=events.dl.p[i] + sum(!is.na(nb_incomplete.placebo[[k]]))
    events.dl.p[i]=sum(nb_incomplete.placebo[[1]][nb_incomplete.placebo[[2]]>0])
    mr.dl.p[i]=events.dl.p[i]/sum(nb_incomplete.placebo[[2]][nb_incomplete.placebo[[2]]>0])
    
    #for (k in 1:n) events.dl.t[i]=events.dl.t[i] + sum(!is.na(nb_incomplete.treatment[[k]]))
    events.dl.t[i]=sum(nb_incomplete.treatment[[1]][nb_incomplete.treatment[[2]]>0])
    mr.dl.t[i]=events.dl.t[i]/sum(nb_incomplete.treatment[[2]][nb_incomplete.treatment[[2]]>0])
    
    #Complete data    
    Yco <- c(cc.p, cc.t)
    alldata=rbind(alldata, data.frame(sim=i, ID=id, TRTn=X, TRT=Xchar, COUNT=Y, cc=Yco, STUDY_TIME=censor.time))
    
    co_poisson.summary <- summary(glm(Yco~X, family=poisson()))
    #Roberts formula for s.e.
    co_poisson.se[i] <- co_poisson.summary$coeffi[2,2]*exp(co_poisson.summary$coeffi[2,1])
    co_poisson.se2[i] <- co_poisson.summary$coeffi[2,2]
    co_poisson.rratio[i] <- exp(co_poisson.summary$coeffi[2,1])   #Rate ratio
    co_poisson.pval[i] <- co_poisson.summary$coef[2,4] #p-value
    
    co_quasipoisson.summary <- summary(glm(Yco~X, family=quasipoisson()))
    #Roberts formula for s.e.
    co_quasipoisson.se[i] <- co_quasipoisson.summary$coeffi[2,2]*exp(co_quasipoisson.summary$coeffi[2,1])
    co_quasipoisson.se2[i] <- co_quasipoisson.summary$coeffi[2,2]
    co_quasipoisson.rratio[i] <- exp(co_quasipoisson.summary$coeffi[2,1])   
    co_quasipoisson.pval[i] <- co_quasipoisson.summary$coef[2,4] 
    
    co_nbd.summary <- summary(glm.nb(Yco~X))
    co_nbd.warn[i] <- (!is.null(co_nbd.summary$th.warn))
    #Roberts formula for s.e.
    co_nbd.se[i] <- co_nbd.summary$coeffi[2,2]*exp(co_nbd.summary$coeffi[2,1])
    co_nbd.se2[i] <- co_nbd.summary$coeffi[2,2]
    co_nbd.rratio[i] <- exp(co_nbd.summary$coeffi[2,1])   
    co_nbd.pval[i] <- co_nbd.summary$coef[2,4] 
    co_nbd.theta[i] <- co_nbd.summary$theta
    
    # no imputation (direct likelihood approach)
    dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    #Roberts formula for s.e.
    dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]*exp(dl_poisson.summary$coeffi[2,1])
    dl_poisson.se2[i] <- dl_poisson.summary$coeffi[2,2]
    dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
    dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
    
    dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    
    #Roberts formula for s.e.
    dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]*exp(dl_quasipoisson.summary$coeffi[2,1])
    dl_quasipoisson.se2[i] <- dl_quasipoisson.summary$coeffi[2,2]
    dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
    dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
    
    dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
    #Roberts formula for s.e.
    dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]*exp(dl_nbd.summary$coeffi[2,1])
    dl_nbd.se2[i] <- dl_nbd.summary$coeffi[2,2]
    dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
    dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
    dl_nbd.theta[i] <- dl_nbd.summary$theta
    
    # jump to reference
    
    j2r <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    j2r_poisson.se[i] <- (j2r$mi_poisson.se2)*exp(j2r$mi_poisson.trt2)
    j2r_poisson.se2[i] <- j2r$mi_poisson.se2
    j2r_poisson.rratio[i] <- exp(j2r$mi_poisson.trt2)
    j2r_poisson.df[i] <- j2r$mi_poisson.df2
    j2r_poisson.adf[i] <- j2r$mi_poisson.adf2
    j2r_poisson.pval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se2[i]),df=j2r_poisson.df[i]))
    j2r_poisson.apval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se2[i]),df=j2r_poisson.adf[i]))
    
    j2r_quasipoisson.se[i] <- (j2r$mi_quasipoisson.se2)*exp(j2r$mi_quasipoisson.trt2)
    j2r_quasipoisson.se2[i] <- j2r$mi_quasipoisson.se2
    j2r_quasipoisson.rratio[i] <- exp(j2r$mi_quasipoisson.trt2)
    j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df2
    j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf2
    j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se2[i]),df=j2r_quasipoisson.df[i]))
    j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se2[i]),df=j2r_quasipoisson.adf[i]))
    
    j2r_nbd.se[i] <- (j2r$mi_nbd.se2)*exp(j2r$mi_nbd.trt2)
    j2r_nbd.se2[i] <- j2r$mi_nbd.se2
    j2r_nbd.rratio[i] <- exp(j2r$mi_nbd.trt2)
    j2r_nbd.df[i] <- j2r$mi_nbd.df2
    j2r_nbd.adf[i] <- j2r$mi_nbd.adf2
    j2r_nbd.pval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se2[i]),df=j2r_nbd.df[i]))
    j2r_nbd.apval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se2[i]),df=j2r_nbd.adf[i]))
    
    mr.j2r.p[i]=j2r$mr.p
    events.j2r.p[i]=j2r$ev.p
    
    mr.j2r.t[i]=j2r$mr.t
    events.j2r.t[i]=j2r$ev.t
    
    #if (!is.null(j2r$outliers)) outliers.j2r=rbind(outliers.j2r, cbind(j2r$outliers,i))
    #print(cbind(j2r$outliers,i))
    
    # efficacy
    
    efy <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)
    efy_poisson.se[i] <- (efy$mi_poisson.se2)*exp(efy$mi_poisson.trt2)
    efy_poisson.se2[i] <- efy$mi_poisson.se2
    efy_poisson.rratio[i] <- exp(efy$mi_poisson.trt2)
    efy_poisson.df[i] <- efy$mi_poisson.df2
    efy_poisson.adf[i] <- efy$mi_poisson.adf2
    efy_poisson.pval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se2[i]),df=efy_poisson.df[i]))
    efy_poisson.apval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se2[i]),df=efy_poisson.adf[i]))
    
    efy_quasipoisson.se[i] <- (efy$mi_quasipoisson.se2)*exp(efy$mi_quasipoisson.trt2)
    efy_quasipoisson.se2[i] <- efy$mi_quasipoisson.se2
    efy_quasipoisson.rratio[i] <- exp(efy$mi_quasipoisson.trt2)
    efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df2
    efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf2
    efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se2[i]),df=efy_quasipoisson.df[i]))
    efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se2[i]),df=efy_quasipoisson.adf[i]))
    
    efy_nbd.se[i] <- (efy$mi_nbd.se2)*exp(efy$mi_nbd.trt2)
    efy_nbd.se2[i] <- efy$mi_nbd.se2
    efy_nbd.rratio[i] <- exp(efy$mi_nbd.trt2)
    efy_nbd.df[i] <- efy$mi_nbd.df2
    efy_nbd.adf[i] <- efy$mi_nbd.adf2
    efy_nbd.pval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se2[i]),df=efy_nbd.df[i]))
    efy_nbd.apval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se2[i]),df=efy_nbd.adf[i]))
    
    mr.ef.p[i]=efy$mr.p
    events.ef.p[i]=efy$ev.p
    mr.ef.t[i]=efy$mr.t
    events.ef.t[i]=efy$ev.t
    
    #if (!is.null(efy$outliers)) outliers.ef=rbind(outliers.ef, cbind(efy$outliers, i))
    
    p.rate[i]=efy$p.rate
    t.rate[i]=efy$t.rate
    gamma.p[i]=efy$gamma.p
    gamma.t[i]=efy$gamma.t
    p.rate2[i]=efy$p.rate2
    t.rate2[i]=efy$t.rate2
    gamma[i]=efy$gamma
    nbd.warn.p[i]=j2r$warn.p
    nbd.warn.t[i]=j2r$warn.t
    #print(i)
    err[i] <- 0
    }, warning = function(w){return(1)})
  }
 
  # excludes cases with warnings
  co_powerpois<-sum(subset(co_poisson.pval,err==0)<alpha)/length(subset(co_poisson.pval,err==0))
  co_powerqpois<-sum(subset(co_quasipoisson.pval,err==0)<alpha)/length(subset(co_quasipoisson.pval,err==0))
  co_powernb<-sum(subset(co_nbd.pval,err==0)<alpha)/length(subset(co_nbd.pval,err==0))

  dl_powerpois<-sum(subset(dl_poisson.pval,err==0)<alpha)/length(subset(dl_poisson.pval,err==0)) 
  dl_powerqpois<-sum(subset(dl_quasipoisson.pval,err==0)<alpha)/length(subset(dl_quasipoisson.pval,err==0)) 
  dl_powernb<-sum(subset(dl_nbd.pval,err==0)<alpha)/length(subset(dl_nbd.pval,err==0)) 
  
  j2r_powerpois<-sum(subset(j2r_poisson.pval,err==0)<alpha)/length(subset(j2r_poisson.pval,err==0))
  j2r_powerqpois<-sum(subset(j2r_quasipoisson.pval,err==0)<alpha)/length(subset(j2r_quasipoisson.pval,err==0)) 
  j2r_powernb<-sum(subset(j2r_nbd.pval,err==0)<alpha)/length(subset(j2r_nbd.pval,err==0)) 
  
  j2r_powerpoisa<-sum(subset(j2r_poisson.apval,err==0)<alpha)/length(subset(j2r_poisson.apval,err==0))
  j2r_powerqpoisa<-sum(subset(j2r_quasipoisson.apval,err==0)<alpha)/length(subset(j2r_quasipoisson.apval,err==0))
  j2r_powernba<-sum(subset(j2r_nbd.apval,err==0)<alpha)/length(subset(j2r_nbd.apval,err==0))
  
  efy_powerpois<-sum(subset(efy_poisson.pval,err==0)<alpha)/length(subset(efy_poisson.pval,err==0))
  efy_powerqpois<-sum(subset(efy_quasipoisson.pval,err==0)<alpha)/length(subset(efy_quasipoisson.pval,err==0))
  efy_powernb<-sum(subset(efy_nbd.pval,err==0)<alpha)/length(subset(efy_nbd.pval,err==0))
  
  efy_powerpoisa<-sum(subset(efy_poisson.apval,err==0)<alpha)/length(subset(efy_poisson.apval,err==0))
  efy_powerqpoisa<-sum(subset(efy_quasipoisson.apval,err==0)<alpha)/length(subset(efy_quasipoisson.apval,err==0))
  efy_powernba<-sum(subset(efy_nbd.apval,err==0)<alpha)/length(subset(efy_nbd.apval,err==0))
  
  power <- as.data.frame(cbind(co_powerpois,co_powerqpois,co_powernb,
                               dl_powerpois,dl_powerqpois,dl_powernb,
                               j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba))
  
  rratio <- as.data.frame(cbind(co_poisson.rratio,co_quasipoisson.rratio,co_nbd.rratio,
                                dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,
                                j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,
                                efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio))
  
  se <- as.data.frame(cbind(co_poisson.se,co_quasipoisson.se,co_nbd.se,
                            dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,
                            j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,
                            efy_poisson.se,efy_quasipoisson.se,efy_nbd.se,
                            co_poisson.se2,co_quasipoisson.se2,co_nbd.se2,
                            dl_poisson.se2,dl_quasipoisson.se2,dl_nbd.se2,
                            j2r_poisson.se2,j2r_quasipoisson.se2,j2r_nbd.se2,
                            efy_poisson.se2,efy_quasipoisson.se2,efy_nbd.se2))
  
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,
                               efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf))
  
  dl_summnbtheta<-summary(dl_nbd.theta)
  
  dl_nbthetawarn <- sum(dl_nbd.warn)
  nbthetawarn.p <- sum(nbd.warn.p)
  nbthetawarn.t <- sum(nbd.warn.t)
  warnalltype <- sum(err==1)
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
  
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  
  write.csv(alldata,datname)
  
  #paste("Data saved as",datname)
  
  result <- list(err=err,power=power,dl_theta=dl_summnbtheta,rratio=rratio,se=se,mi_df=mi_df,dl_nbthetawarn=dl_nbthetawarn,dropout=dropout,censortime=censortime,ERgt75p=ERgt75p,ERle75p=ERle75p, py.p=py.p, py.t=py.t, 
                 mr.dl.p=mr.dl.p, mr.ef.p=mr.ef.p, mr.co.p=mr.co.p, mr.j2r.p=mr.j2r.p, events.dl.p=events.dl.p, events.ef.p=events.ef.p, events.j2r.p=events.j2r.p, events.co.p=events.co.p,
                 mr.dl.t=mr.dl.t, mr.ef.t=mr.ef.t, mr.co.t=mr.co.t, mr.j2r.t=mr.j2r.t, events.dl.t=events.dl.t, events.ef.t=events.ef.t, events.j2r.t=events.j2r.t, events.co.t=events.co.t, 
                 p.rate=p.rate, t.rate=t.rate, p.rate2=p.rate2, t.rate2=t.rate2, gamma.p=gamma.p, gamma.t=gamma.t, gamma=gamma, nbthetawarn.p=nbthetawarn.p, nbthetawarn.t=nbthetawarn.t, warnalltype=warnalltype,outliers.ef=outliers.ef, outliers.j2r=outliers.j2r)
  return(result)
}



sim_mar2_mg<-function(n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,drop_rate_variance,N.rep,N.mi, alpha, datname){
  
  rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  
  py.p <- rep(NA,N.rep)
  py.t <- rep(NA,N.rep)
  
  events.dl.p <- rep(NA,N.rep)
  mr.dl.p <- rep(NA,N.rep)
  mr.ef.p <- rep(NA,N.rep)
  events.ef.p <- rep(NA,N.rep)
  mr.co.p <- rep(NA,N.rep)
  events.co.p <- rep(NA,N.rep)
  mr.j2r.p <- rep(NA,N.rep)
  events.j2r.p <- rep(NA,N.rep)
  
  
  events.dl.t <- rep(NA,N.rep)
  mr.dl.t <- rep(NA,N.rep)
  mr.ef.t <- rep(NA,N.rep)
  events.ef.t <- rep(NA,N.rep)
  mr.co.t <- rep(NA,N.rep)
  events.co.t <- rep(NA,N.rep)
  mr.j2r.t <- rep(NA,N.rep)
  events.j2r.t <- rep(NA,N.rep)
  
  censor.time.p <- rep(NA,N.rep)
  censor.time.t <- rep(NA,N.rep)
  t.rate <- rep(NA,N.rep)
  p.rate <- rep(NA,N.rep)
  
  ERgt75p.p <- rep(NA,N.rep)
  ERgt75p.t <- rep(NA,N.rep)
  ERle75p.p <- rep(NA,N.rep)
  ERle75p.t <- rep(NA,N.rep)
  
  # no imputation (direct likelihood approach)
  dl_poisson.se <- rep(NA,N.rep)
  dl_poisson.se2 <- rep(NA,N.rep)
  dl_poisson.rratio <- rep(NA,N.rep)
  dl_poisson.pval <- rep(NA,N.rep)
  
  dl_quasipoisson.se <- rep(NA,N.rep)
  dl_quasipoisson.se2 <- rep(NA,N.rep)
  dl_quasipoisson.rratio <- rep(NA,N.rep)
  dl_quasipoisson.pval <- rep(NA,N.rep)
  
  dl_nbd.se <- rep(NA,N.rep)
  dl_nbd.se2 <- rep(NA,N.rep)
  dl_nbd.rratio <- rep(NA,N.rep)
  dl_nbd.pval <- rep(NA,N.rep)
  dl_nbd.theta <- rep(NA,N.rep)
  dl_nbd.warn <- rep(NA,N.rep)
  nbd.warn.p <- rep(NA,N.rep)
  nbd.warn.t <- rep(NA,N.rep)
  
  # jump to reference 
  j2r_poisson.se <- rep(NA,N.rep)
  j2r_poisson.se2 <- rep(NA,N.rep)
  j2r_poisson.rratio <- rep(NA,N.rep)
  j2r_poisson.pval <- rep(NA,N.rep)
  j2r_poisson.apval <- rep(NA,N.rep)
  j2r_poisson.df <- rep(NA,N.rep)
  j2r_poisson.adf <- rep(NA,N.rep) 
  
  j2r_quasipoisson.se <- rep(NA,N.rep)
  j2r_quasipoisson.se2 <- rep(NA,N.rep)
  j2r_quasipoisson.rratio <- rep(NA,N.rep)
  j2r_quasipoisson.pval <- rep(NA,N.rep)
  j2r_quasipoisson.apval <- rep(NA,N.rep)
  j2r_quasipoisson.df <- rep(NA,N.rep)
  j2r_quasipoisson.adf <- rep(NA,N.rep)
  
  j2r_nbd.se <- rep(NA,N.rep)
  j2r_nbd.se2 <- rep(NA,N.rep)
  j2r_nbd.rratio <- rep(NA,N.rep)
  j2r_nbd.pval <- rep(NA,N.rep)
  j2r_nbd.apval <- rep(NA,N.rep)
  j2r_nbd.df <- rep(NA,N.rep)
  j2r_nbd.adf <- rep(NA,N.rep)
  outliers.j2r=c()
  
  # efficacy estimand
  efy_poisson.se <- rep(NA,N.rep)
  efy_poisson.se2 <- rep(NA,N.rep)
  efy_poisson.rratio <- rep(NA,N.rep)
  efy_poisson.pval <- rep(NA,N.rep)
  efy_poisson.apval <- rep(NA,N.rep)
  efy_poisson.df <- rep(NA,N.rep)
  efy_poisson.adf <- rep(NA,N.rep)
  
  efy_quasipoisson.se <- rep(NA,N.rep)
  efy_quasipoisson.se2 <- rep(NA,N.rep)
  efy_quasipoisson.rratio <- rep(NA,N.rep)
  efy_quasipoisson.pval <- rep(NA,N.rep)
  efy_quasipoisson.apval <- rep(NA,N.rep)
  efy_quasipoisson.df <- rep(NA,N.rep)
  efy_quasipoisson.adf <- rep(NA,N.rep)
  
  efy_nbd.se <- rep(NA,N.rep)
  efy_nbd.se2 <- rep(NA,N.rep)
  efy_nbd.rratio <- rep(NA,N.rep)
  efy_nbd.pval <- rep(NA,N.rep)
  efy_nbd.apval <- rep(NA,N.rep)
  efy_nbd.df <- rep(NA,N.rep)
  efy_nbd.adf <- rep(NA,N.rep)
  outliers.ef=c()
  alldata=data.frame()
  
  for (i in 1:N.rep){
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- dropout_time4(T,n,nb_complete.placebo,rates,drop_rate_variance)
    dropout_time.treatment <- dropout_time4(T,n,nb_complete.treatment,rates,drop_rate_variance)
    
    dropout.p[i] <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t[i] <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    cc.p=rep(0,n)
    for (k in 1:n) cc.p[k]=sum(!is.na(nb_complete.placebo[[k]]))
    cc.t=rep(0,n)
    for (k in 1:n) cc.t[k]=sum(!is.na(nb_complete.treatment[[k]]))
    
    alldata=rbind(alldata, cbind(makeDataFrame(nb_incomplete.placebo,nb_incomplete.treatment),cc=c(cc.p,cc.t),sim=i))
    
    # mean ER of patients who stayed for at least 75% of total follow-up period vs early dropouts
    ERgt75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]>0.75*T))
    ERgt75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]>0.75*T))
    ERle75p.p[i] <- mean(subset(nb_incomplete.placebo[[1]]/nb_incomplete.placebo[[2]],nb_incomplete.placebo[[2]]<=0.75*T))
    ERle75p.t[i] <- mean(subset(nb_incomplete.treatment[[1]]/nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<=0.75*T))
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p[i] <- mean(subset(nb_incomplete.placebo[[2]],  nb_incomplete.placebo[[2]]<T))
    censor.time.t[i] <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    #MG add
    py.p[i]=sum(nb_incomplete.placebo[[2]])
    py.t[i]=sum(nb_incomplete.treatment[[2]])
    
    events.co.p[i]=0
    for (k in 1:n) events.co.p[i]=events.co.p[i] + sum(!is.na(nb_complete.placebo[[k]]))                 
    mr.co.p[i]=events.co.p[i]/(n)
    
    events.co.t[i]=0
    for (k in 1:n) events.co.t[i]=events.co.t[i] + sum(!is.na(nb_complete.treatment[[k]]))                
    mr.co.t[i]=events.co.t[i]/(n)
    
    #for (k in 1:n) events.dl.p[i]=events.dl.p[i] + sum(!is.na(nb_incomplete.placebo[[k]]))
    events.dl.p[i]=sum(nb_incomplete.placebo[[1]][nb_incomplete.placebo[[2]]>0])
    mr.dl.p[i]=events.dl.p[i]/sum(nb_incomplete.placebo[[2]][nb_incomplete.placebo[[2]]>0])
                                                              
    #for (k in 1:n) events.dl.t[i]=events.dl.t[i] + sum(!is.na(nb_incomplete.treatment[[k]]))
    events.dl.t[i]=sum(nb_incomplete.treatment[[1]][nb_incomplete.treatment[[2]]>0])
    mr.dl.t[i]=events.dl.t[i]/sum(nb_incomplete.treatment[[2]][nb_incomplete.treatment[[2]]>0])
    
    # no imputation (direct likelihood approach)
    
    dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    #Roberts formula for s.e.
    dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]*exp(dl_poisson.summary$coeffi[2,1])
    dl_poisson.se2[i] <- dl_poisson.summary$coeffi[2,2]
    dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
    dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
    
    dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    
    #Roberts formula for s.e.
    dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]*exp(dl_quasipoisson.summary$coeffi[2,1])
    dl_quasipoisson.se2[i] <- dl_quasipoisson.summary$coeffi[2,2]
    dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
    dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
    
    dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
    #Roberts formula for s.e.
    dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]*exp(dl_nbd.summary$coeffi[2,1])
    dl_nbd.se2[i] <- dl_nbd.summary$coeffi[2,2]
    dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
    dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
    dl_nbd.theta[i] <- dl_nbd.summary$theta
    
    # jump to reference
    
    j2r <- cim_mg(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    j2r_poisson.se[i] <- j2r$mi_poisson.se
    j2r_poisson.se2[i] <- j2r$mi_poisson.se2
    j2r_poisson.rratio[i] <- j2r$mi_poisson.trt
    j2r_poisson.df[i] <- j2r$mi_poisson.df
    j2r_poisson.adf[i] <- j2r$mi_poisson.adf
    j2r_poisson.pval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se2[i]),df=j2r_poisson.df[i]))
    j2r_poisson.apval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se2[i]),df=j2r_poisson.adf[i]))
    
    j2r_quasipoisson.se[i] <- j2r$mi_quasipoisson.se
    j2r_quasipoisson.se2[i] <- j2r$mi_quasipoisson.se2
    j2r_quasipoisson.rratio[i] <- j2r$mi_quasipoisson.trt
    j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df
    j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf
    j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se2[i]),df=j2r_quasipoisson.df[i]))
    j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se2[i]),df=j2r_quasipoisson.adf[i]))
    
    j2r_nbd.se[i] <- j2r$mi_nbd.se
    j2r_nbd.se2[i] <- j2r$mi_nbd.se2
    j2r_nbd.rratio[i] <- j2r$mi_nbd.trt
    j2r_nbd.df[i] <- j2r$mi_nbd.df
    j2r_nbd.adf[i] <- j2r$mi_nbd.adf
    j2r_nbd.pval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se2[i]),df=j2r_nbd.df[i]))
    j2r_nbd.apval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se2[i]),df=j2r_nbd.adf[i]))
    
    mr.j2r.p[i]=j2r$mr.p
    events.j2r.p[i]=j2r$ev.p
    
    mr.j2r.t[i]=j2r$mr.t
    events.j2r.t[i]=j2r$ev.t
    
    if (!is.null(j2r$outliers)) outliers.j2r=rbind(outliers.j2r, cbind(j2r$outliers,i))
    #print(cbind(j2r$outliers,i))
    
    # efficacy
    
    efy <- cim_mg(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)
    efy_poisson.se[i] <- efy$mi_poisson.se
    efy_poisson.se2[i] <- efy$mi_poisson.se2
    efy_poisson.rratio[i] <- efy$mi_poisson.trt
    efy_poisson.df[i] <- efy$mi_poisson.df
    efy_poisson.adf[i] <- efy$mi_poisson.adf
    efy_poisson.pval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se2[i]),df=efy_poisson.df[i]))
    efy_poisson.apval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se2[i]),df=efy_poisson.adf[i]))
    
    efy_quasipoisson.se[i] <- efy$mi_quasipoisson.se
    efy_quasipoisson.se2[i] <- efy$mi_quasipoisson.se2
    efy_quasipoisson.rratio[i] <- efy$mi_quasipoisson.trt
    efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df
    efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf
    efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se2[i]),df=efy_quasipoisson.df[i]))
    efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se2[i]),df=efy_quasipoisson.adf[i]))
    
    efy_nbd.se[i] <- efy$mi_nbd.se
    efy_nbd.se2[i] <- efy$mi_nbd.se2
    efy_nbd.rratio[i] <- efy$mi_nbd.trt
    efy_nbd.df[i] <- efy$mi_nbd.df
    efy_nbd.adf[i] <- efy$mi_nbd.adf
    efy_nbd.pval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se2[i]),df=efy_nbd.df[i]))
    efy_nbd.apval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se2[i]),df=efy_nbd.adf[i]))
    
    mr.ef.p[i]=efy$mr.p
    events.ef.p[i]=efy$ev.p
    mr.ef.t[i]=efy$mr.t
    events.ef.t[i]=efy$ev.t
    
    if (!is.null(efy$outliers)) outliers.ef=rbind(outliers.ef, cbind(efy$outliers, i))
      
    p.rate[i]=efy$p.rate
    t.rate[i]=efy$t.rate
    nbd.warn.p[i]=j2r$warn.p
    nbd.warn.t[i]=j2r$warn.t
    print(i)
  }
  
  dl_powerpois<-sum(dl_poisson.pval<alpha)/N.rep 
  dl_powerqpois<-sum(dl_quasipoisson.pval<alpha)/N.rep 
  dl_powernb<-sum(dl_nbd.pval<alpha)/N.rep 
  
  j2r_powerpois<-sum(j2r_poisson.pval<alpha)/N.rep 
  j2r_powerqpois<-sum(j2r_quasipoisson.pval<alpha)/N.rep 
  j2r_powernb<-sum(j2r_nbd.pval<alpha)/N.rep 
  
  j2r_powerpoisa<-sum(j2r_poisson.apval<alpha)/N.rep 
  j2r_powerqpoisa<-sum(j2r_quasipoisson.apval<alpha)/N.rep 
  j2r_powernba<-sum(j2r_nbd.apval<alpha)/N.rep 
  
  efy_powerpois<-sum(efy_poisson.pval<alpha)/N.rep 
  efy_powerqpois<-sum(efy_quasipoisson.pval<alpha)/N.rep 
  efy_powernb<-sum(efy_nbd.pval<alpha)/N.rep 
  
  efy_powerpoisa<-sum(efy_poisson.apval<alpha)/N.rep 
  efy_powerqpoisa<-sum(efy_quasipoisson.apval<alpha)/N.rep 
  efy_powernba<-sum(efy_nbd.apval<alpha)/N.rep 
  
  power <- as.data.frame(cbind(dl_powerpois,dl_powerqpois,dl_powernb,j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba))
  
  rratio <- as.data.frame(cbind(dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio))
  se <- as.data.frame(cbind(dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,efy_poisson.se,efy_quasipoisson.se,efy_nbd.se,
                            dl_poisson.se2,dl_quasipoisson.se2,dl_nbd.se2,j2r_poisson.se2,j2r_quasipoisson.se2,j2r_nbd.se2,efy_poisson.se2,efy_quasipoisson.se2,efy_nbd.se2))
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf))
  
  dl_summnbtheta<-summary(dl_nbd.theta)
  
  dl_nbthetawarn <- sum(dl_nbd.warn)
  nbthetawarn.p <- sum(nbd.warn.p)
  nbthetawarn.t <- sum(nbd.warn.t)
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
 
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  write.csv(alldata,datname)
  
  paste("Data saved as",datname)
  
  result <- list(power=power,dl_theta=dl_summnbtheta,rratio=rratio,se=se,mi_df=mi_df,dl_nbthetawarn=dl_nbthetawarn,dropout=dropout,censortime=censortime,ERgt75p=ERgt75p,ERle75p=ERle75p, py.p=py.p, py.t=py.t, 
                 mr.dl.p=mr.dl.p, mr.ef.p=mr.ef.p, mr.co.p=mr.co.p, mr.j2r.p=mr.j2r.p, events.dl.p=events.dl.p, events.ef.p=events.ef.p, events.j2r.p=events.j2r.p, events.co.p=events.co.p,
                 mr.dl.t=mr.dl.t, mr.ef.t=mr.ef.t, mr.co.t=mr.co.t, mr.j2r.t=mr.j2r.t, events.dl.t=events.dl.t, events.ef.t=events.ef.t, events.j2r.t=events.j2r.t, events.co.t=events.co.t, 
                 p.rate=p.rate, t.rate=t.rate, outliers.ef=outliers.ef, outliers.j2r=outliers.j2r, nbthetawarn.p=nbthetawarn.p, nbthetawarn.t=nbthetawarn.t)
  return(result)
}




gen1mar=function(n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,drop_rate_variance,sasname){
  
  rate.treatment <- (1-rr)*rate.placebo
  
  nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
  nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
  
  dropout_time.placebo <- dropout_time4(T,n,nb_complete.placebo,rates,drop_rate_variance)
  dropout_time.treatment <- dropout_time4(T,n,nb_complete.treatment,rates,drop_rate_variance)
  
  dropout.p <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
  dropout.t <- length(dropout_time.treatment[dropout_time.treatment < T])/n
  
  nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
  nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
  
  Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
  X <- c(rep(0,n),rep(1,n))
  censor.time.p <- mean(subset(nb_incomplete.placebo[[2]],  nb_incomplete.placebo[[2]]<T))
  censor.time.t <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
  censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
  
  id=1:n
  Xchar=rep("",length(X)) 
  Xchar[X==0]="Placebo"
  Xchar[X==1]="Treatment"
  
  data=data.frame(ID=id,TRTn=X,TRT=Xchar,COUNT=Y,STUDY_TIME=censor.time)
  write.csv(data,sasname)
  return(data)
}


genNmar=function(N,n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,drop_rate_variance){
  data=data.frame()
  for (i in 1:N) {
    
    rate.treatment <- (1-rr)*rate.placebo
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- dropout_time4(T,n,nb_complete.placebo,rates,drop_rate_variance)
    dropout_time.treatment <- dropout_time4(T,n,nb_complete.treatment,rates,drop_rate_variance)
    
    dropout.p <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p <- mean(subset(nb_incomplete.placebo[[2]],  nb_incomplete.placebo[[2]]<T))
    censor.time.t <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    id=1:n
    Xchar=rep("",length(X)) 
    Xchar[X==0]="Placebo"
    Xchar[X==1]="Treatment"
    
    data=rbind(data, data.frame(sim=i, ID=id, TRTn=X, TRT=Xchar, COUNT=Y, STUDY_TIME=censor.time))
  }
  return(data)
}

genNmar.mi=function(N,n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,drop_rate_variance, N.mi){
  data=data.frame()
  for (i in 1:N) {
    
    rate.treatment <- (1-rr)*rate.placebo
    
    nb_complete.placebo <- nb_complete(T,n,rate.placebo,dispersion.placebo)
    nb_complete.treatment <- nb_complete(T,n,rate.treatment,dispersion.treatment)
    
    dropout_time.placebo <- dropout_time4(T,n,nb_complete.placebo,rates,drop_rate_variance)
    dropout_time.treatment <- dropout_time4(T,n,nb_complete.treatment,rates,drop_rate_variance)
    
    dropout.p <- length(dropout_time.placebo[dropout_time.placebo < T])/n  
    dropout.t <- length(dropout_time.treatment[dropout_time.treatment < T])/n
    
    nb_incomplete.placebo <- nb_incomplete(n,nb_complete.placebo,dropout_time.placebo)
    nb_incomplete.treatment <- nb_incomplete(n,nb_complete.treatment,dropout_time.treatment)
    
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time.p <- mean(subset(nb_incomplete.placebo[[2]],  nb_incomplete.placebo[[2]]<T))
    censor.time.t <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    
    id=1:n
    Xchar=rep("",length(X)) 
    Xchar[X==0]="Placebo"
    Xchar[X==1]="Treatment"
    
    Imp.j2r = cim_imponly(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    Imp.efy = cim_imponly(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)    
    data=rbind(data, data.frame(sim=i, ID=id, TRTn=X, TRT=Xchar, COUNT=Y, STUDY_TIME=censor.time, Imp.id=Imp.j2r$imp, Imp.j2r=Imp.j2r$Imputed, Imp.efy=Imp.efy$Imputed))
  }
  return(data)
}

csumev = function(data, N.sim, N.imp, trt){
  mj2r=ddply(data, c("sim","ID","TRT"), function(df) mean(df$Imp.j2r))
  mefy=ddply(data, c("sim","ID","TRT"), function(df) mean(df$Imp.efy))
  outcum=data.frame()
  for (i in 1:N.sim) {
    for (trti in trt) {
        dati=data[data$TRT==trti & data$sim==i & data$Imp.id==1,]
        mj2ri=mj2r[mj2r$TRT==trti & mj2r$sim==i,]
        mefyi=mefy[mefy$TRT==trti & mefy$sim==i,]
        cc=rbind(data.frame(Method="Incomplete", id=dati$ID, time=dati$STUDY_TIME, count=(dati$COUNT), ccount=cumsum(dati$COUNT),  trt=trti,sim=i), 
                 data.frame(Method="J2R"       , id=dati$ID, time=dati$STUDY_TIME, count=(mj2ri$V1), ccount=cumsum(mj2ri$V1),trt=trti,sim=i),
                 data.frame(Method="MAR"       , id=dati$ID, time=dati$STUDY_TIME, count=(mefyi$V1), ccount=cumsum(mefyi$V1),trt=trti,sim=i))
        outcum=rbind(outcum, cc)
    }
  }
  return(outcum)
}



# MAR with variation in drop-out rate
# Read data from file
sim_mar2_rd<-function(n,N.rep,N.mi, alpha,T, datname){
  
  
  #rate.treatment <- (1-rr)*rate.placebo
  dropout.p <- rep(NA,N.rep)
  dropout.t <- rep(NA,N.rep)
  
  py.p <- rep(NA,N.rep)
  py.t <- rep(NA,N.rep)
  
  events.dl.p <- rep(NA,N.rep)
  mr.dl.p <- rep(NA,N.rep)
  mr.ef.p <- rep(NA,N.rep)
  events.ef.p <- rep(NA,N.rep)
  mr.co.p <- rep(NA,N.rep)
  events.co.p <- rep(NA,N.rep)
  mr.j2r.p <- rep(NA,N.rep)
  events.j2r.p <- rep(NA,N.rep)
  
  
  events.dl.t <- rep(NA,N.rep)
  mr.dl.t <- rep(NA,N.rep)
  mr.ef.t <- rep(NA,N.rep)
  events.ef.t <- rep(NA,N.rep)
  mr.co.t <- rep(NA,N.rep)
  events.co.t <- rep(NA,N.rep)
  mr.j2r.t <- rep(NA,N.rep)
  events.j2r.t <- rep(NA,N.rep)
  
  censor.time.p <- rep(NA,N.rep)
  censor.time.t <- rep(NA,N.rep)
  t.rate <- rep(NA,N.rep)
  p.rate <- rep(NA,N.rep)
  gamma.t <- rep(NA,N.rep)
  gamma.p <- rep(NA,N.rep)
  
  ERgt75p.p <- rep(NA,N.rep)
  ERgt75p.t <- rep(NA,N.rep)
  ERle75p.p <- rep(NA,N.rep)
  ERle75p.t <- rep(NA,N.rep)
  
  #Complete data
  co_poisson.se <- rep(NA,N.rep)
  co_poisson.se2 <- rep(NA,N.rep)
  co_poisson.rratio <- rep(NA,N.rep)
  co_poisson.pval <- rep(NA,N.rep)
  
  co_quasipoisson.se <- rep(NA,N.rep)
  co_quasipoisson.se2 <- rep(NA,N.rep)
  co_quasipoisson.rratio <- rep(NA,N.rep)
  co_quasipoisson.pval <- rep(NA,N.rep)
  
  co_nbd.se <- rep(NA,N.rep)
  co_nbd.se2 <- rep(NA,N.rep)
  co_nbd.rratio <- rep(NA,N.rep)
  co_nbd.pval <- rep(NA,N.rep)
  co_nbd.theta <- rep(NA,N.rep)
  co_nbd.warn <- rep(NA,N.rep)
  
  # no imputation (direct likelihood approach)
  dl_poisson.se <- rep(NA,N.rep)
  dl_poisson.se2 <- rep(NA,N.rep)
  dl_poisson.rratio <- rep(NA,N.rep)
  dl_poisson.pval <- rep(NA,N.rep)
  
  dl_quasipoisson.se <- rep(NA,N.rep)
  dl_quasipoisson.se2 <- rep(NA,N.rep)
  dl_quasipoisson.rratio <- rep(NA,N.rep)
  dl_quasipoisson.pval <- rep(NA,N.rep)
  
  dl_nbd.se <- rep(NA,N.rep)
  dl_nbd.se2 <- rep(NA,N.rep)
  dl_nbd.rratio <- rep(NA,N.rep)
  dl_nbd.pval <- rep(NA,N.rep)
  dl_nbd.theta <- rep(NA,N.rep)
  dl_nbd.warn <- rep(NA,N.rep)
  nbd.warn.p <- rep(NA,N.rep)
  nbd.warn.t <- rep(NA,N.rep)
  
  # jump to reference 
  j2r_poisson.se <- rep(NA,N.rep)
  j2r_poisson.se2 <- rep(NA,N.rep)
  j2r_poisson.rratio <- rep(NA,N.rep)
  j2r_poisson.pval <- rep(NA,N.rep)
  j2r_poisson.apval <- rep(NA,N.rep)
  j2r_poisson.df <- rep(NA,N.rep)
  j2r_poisson.adf <- rep(NA,N.rep) 
  
  j2r_quasipoisson.se <- rep(NA,N.rep)
  j2r_quasipoisson.se2 <- rep(NA,N.rep)
  j2r_quasipoisson.rratio <- rep(NA,N.rep)
  j2r_quasipoisson.pval <- rep(NA,N.rep)
  j2r_quasipoisson.apval <- rep(NA,N.rep)
  j2r_quasipoisson.df <- rep(NA,N.rep)
  j2r_quasipoisson.adf <- rep(NA,N.rep)
  
  j2r_nbd.se <- rep(NA,N.rep)
  j2r_nbd.se2 <- rep(NA,N.rep)
  j2r_nbd.rratio <- rep(NA,N.rep)
  j2r_nbd.pval <- rep(NA,N.rep)
  j2r_nbd.apval <- rep(NA,N.rep)
  j2r_nbd.df <- rep(NA,N.rep)
  j2r_nbd.adf <- rep(NA,N.rep)
  outliers.j2r=c()
  
  # efficacy estimand
  efy_poisson.se <- rep(NA,N.rep)
  efy_poisson.se2 <- rep(NA,N.rep)
  efy_poisson.rratio <- rep(NA,N.rep)
  efy_poisson.pval <- rep(NA,N.rep)
  efy_poisson.apval <- rep(NA,N.rep)
  efy_poisson.df <- rep(NA,N.rep)
  efy_poisson.adf <- rep(NA,N.rep)
  
  efy_quasipoisson.se <- rep(NA,N.rep)
  efy_quasipoisson.se2 <- rep(NA,N.rep)
  efy_quasipoisson.rratio <- rep(NA,N.rep)
  efy_quasipoisson.pval <- rep(NA,N.rep)
  efy_quasipoisson.apval <- rep(NA,N.rep)
  efy_quasipoisson.df <- rep(NA,N.rep)
  efy_quasipoisson.adf <- rep(NA,N.rep)
  
  efy_nbd.se <- rep(NA,N.rep)
  efy_nbd.se2 <- rep(NA,N.rep)
  efy_nbd.rratio <- rep(NA,N.rep)
  efy_nbd.pval <- rep(NA,N.rep)
  efy_nbd.apval <- rep(NA,N.rep)
  efy_nbd.df <- rep(NA,N.rep)
  efy_nbd.adf <- rep(NA,N.rep)
  outliers.ef=c()
  alldata=data.frame()
  
  err <- rep(NA,N.rep)
  
  #Get data
  indata=read.csv(datname)

  for (i in 1:N.rep){
     if ((i/10 - round(i/10)) == 0) print(i)
    err[i] <-tryCatch({
    indatai=indata[indata$sim==i,]
    
    cc.p=indatai$cc[indatai$TRT=="Placebo"]
    cc.t=indatai$cc[indatai$TRT=="Active"]
 
    Y <- indatai$COUNT
    X <- indatai$TRTn
    censor.time.p <- indatai$STUDY_TIME[indatai$TRT=="Placebo"]
    censor.time.t <- indatai$STUDY_TIME[indatai$TRT=="Active"]
    censor.time <- indatai$STUDY_TIME

    dropout.p[i] <- length(censor.time.p[censor.time.p < T])/n  
    dropout.t[i] <- length(censor.time.t[censor.time.t < T])/n
    
    #MG add
    py.p[i]=sum(indatai$STUDY_TIME[indatai$TRT=="Placebo"])
    py.t[i]=sum(indatai$STUDY_TIME[indatai$TRT=="Active"])
    
    events.co.p[i]=sum(cc.p)            
    mr.co.p[i]=events.co.p[i]/n
    
    events.co.t[i]=sum(cc.t)             
    mr.co.t[i]=events.co.t[i]/n
    
    events.dl.p[i]=sum(indatai$COUNT[indatai$TRT=="Placebo"])
    mr.dl.p[i]=events.dl.p[i]/py.p[i]
                                                              
    events.dl.t[i]=sum(indatai$COUNT[indatai$TRT=="Active"])
    mr.dl.t[i]=events.dl.t[i]/py.t[i]
    
    Yco <- c(cc.p, cc.t)
    co_poisson.summary <- summary(glm(Yco~X, family=poisson()))
    #Roberts formula for s.e.
    co_poisson.se[i] <- co_poisson.summary$coeffi[2,2]*exp(co_poisson.summary$coeffi[2,1])
    co_poisson.se2[i] <- co_poisson.summary$coeffi[2,2]
    co_poisson.rratio[i] <- exp(co_poisson.summary$coeffi[2,1])   #Rate ratio
    co_poisson.pval[i] <- co_poisson.summary$coef[2,4] #p-value
    
    co_quasipoisson.summary <- summary(glm(Yco~X, family=quasipoisson()))
    #Roberts formula for s.e.
    co_quasipoisson.se[i] <- co_quasipoisson.summary$coeffi[2,2]*exp(co_quasipoisson.summary$coeffi[2,1])
    co_quasipoisson.se2[i] <- co_quasipoisson.summary$coeffi[2,2]
    co_quasipoisson.rratio[i] <- exp(co_quasipoisson.summary$coeffi[2,1])   
    co_quasipoisson.pval[i] <- co_quasipoisson.summary$coef[2,4] 
    
    co_nbd.summary <- summary(glm.nb(Yco~X))
    co_nbd.warn[i] <- (!is.null(co_nbd.summary$th.warn))
    #Roberts formula for s.e.
    co_nbd.se[i] <- co_nbd.summary$coeffi[2,2]*exp(co_nbd.summary$coeffi[2,1])
    co_nbd.se2[i] <- co_nbd.summary$coeffi[2,2]
    co_nbd.rratio[i] <- exp(co_nbd.summary$coeffi[2,1])   
    co_nbd.pval[i] <- co_nbd.summary$coef[2,4] 
    co_nbd.theta[i] <- co_nbd.summary$theta
    
    # no imputation (direct likelihood approach)
    dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    #Roberts formula for s.e.
    dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]*exp(dl_poisson.summary$coeffi[2,1])
    dl_poisson.se2[i] <- dl_poisson.summary$coeffi[2,2]
    dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
    dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
    
    dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    
    #Roberts formula for s.e.
    dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]*exp(dl_quasipoisson.summary$coeffi[2,1])
    dl_quasipoisson.se2[i] <- dl_quasipoisson.summary$coeffi[2,2]
    dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
    dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
    
    dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
    #Roberts formula for s.e.
    dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]*exp(dl_nbd.summary$coeffi[2,1])
    dl_nbd.se2[i] <- dl_nbd.summary$coeffi[2,2]
    dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
    dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
    dl_nbd.theta[i] <- dl_nbd.summary$theta
    
    # jump to reference
    
    j2r <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi)
    j2r_poisson.se[i] <- (j2r$mi_poisson.se2)*exp(j2r$mi_poisson.trt2)
    j2r_poisson.se2[i] <- j2r$mi_poisson.se2
    j2r_poisson.rratio[i] <- exp(j2r$mi_poisson.trt2)
    j2r_poisson.df[i] <- j2r$mi_poisson.df2
    j2r_poisson.adf[i] <- j2r$mi_poisson.adf2
    j2r_poisson.pval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se2[i]),df=j2r_poisson.df[i]))
    j2r_poisson.apval[i] <- 2*(1-pt(abs(log(j2r_poisson.rratio[i])/j2r_poisson.se2[i]),df=j2r_poisson.adf[i]))
    
    j2r_quasipoisson.se[i] <- (j2r$mi_quasipoisson.se2)*exp(j2r$mi_quasipoisson.trt2)
    j2r_quasipoisson.se2[i] <- j2r$mi_quasipoisson.se2
    j2r_quasipoisson.rratio[i] <- exp(j2r$mi_quasipoisson.trt2)
    j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df2
    j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf2
    j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se2[i]),df=j2r_quasipoisson.df[i]))
    j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(j2r_quasipoisson.rratio[i])/j2r_quasipoisson.se2[i]),df=j2r_quasipoisson.adf[i]))
    
    j2r_nbd.se[i] <- (j2r$mi_nbd.se2)*exp(j2r$mi_nbd.trt2)
    j2r_nbd.se2[i] <- j2r$mi_nbd.se2
    j2r_nbd.rratio[i] <- exp(j2r$mi_nbd.trt2)
    j2r_nbd.df[i] <- j2r$mi_nbd.df2
    j2r_nbd.adf[i] <- j2r$mi_nbd.adf2
    j2r_nbd.pval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se2[i]),df=j2r_nbd.df[i]))
    j2r_nbd.apval[i] <- 2*(1-pt(abs(log(j2r_nbd.rratio[i])/j2r_nbd.se2[i]),df=j2r_nbd.adf[i]))
    
    mr.j2r.p[i]=j2r$mr.p
    events.j2r.p[i]=j2r$ev.p
    
    mr.j2r.t[i]=j2r$mr.t
    events.j2r.t[i]=j2r$ev.t
    
    #if (!is.null(j2r$outliers)) outliers.j2r=rbind(outliers.j2r, cbind(j2r$outliers,i))
    #print(cbind(j2r$outliers,i))
    
    # efficacy
    
    efy <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi)
    efy_poisson.se[i] <- (efy$mi_poisson.se2)*exp(efy$mi_poisson.trt2)
    efy_poisson.se2[i] <- efy$mi_poisson.se2
    efy_poisson.rratio[i] <- exp(efy$mi_poisson.trt2)
    efy_poisson.df[i] <- efy$mi_poisson.df2
    efy_poisson.adf[i] <- efy$mi_poisson.adf2
    efy_poisson.pval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se2[i]),df=efy_poisson.df[i]))
    efy_poisson.apval[i] <- 2*(1-pt(abs(log(efy_poisson.rratio[i])/efy_poisson.se2[i]),df=efy_poisson.adf[i]))
    
    efy_quasipoisson.se[i] <- (efy$mi_quasipoisson.se2)*exp(efy$mi_quasipoisson.trt2)
    efy_quasipoisson.se2[i] <- efy$mi_quasipoisson.se2
    efy_quasipoisson.rratio[i] <- exp(efy$mi_quasipoisson.trt2)
    efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df2
    efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf2
    efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se2[i]),df=efy_quasipoisson.df[i]))
    efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(efy_quasipoisson.rratio[i])/efy_quasipoisson.se2[i]),df=efy_quasipoisson.adf[i]))
    
    efy_nbd.se[i] <- (efy$mi_nbd.se2)*exp(efy$mi_nbd.trt2)
    efy_nbd.se2[i] <- efy$mi_nbd.se2
    efy_nbd.rratio[i] <- exp(efy$mi_nbd.trt2)
    efy_nbd.df[i] <- efy$mi_nbd.df2
    efy_nbd.adf[i] <- efy$mi_nbd.adf2
    efy_nbd.pval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se2[i]),df=efy_nbd.df[i]))
    efy_nbd.apval[i] <- 2*(1-pt(abs(log(efy_nbd.rratio[i])/efy_nbd.se2[i]),df=efy_nbd.adf[i]))
    
    mr.ef.p[i]=efy$mr.p
    events.ef.p[i]=efy$ev.p
    mr.ef.t[i]=efy$mr.t
    events.ef.t[i]=efy$ev.t
    
    #if (!is.null(efy$outliers)) outliers.ef=rbind(outliers.ef, cbind(efy$outliers, i))
    
    p.rate[i]=efy$p.rate
    t.rate[i]=efy$t.rate
    gamma.p[i]=efy$gamma.p
    gamma.t[i]=efy$gamma.t
    nbd.warn.p[i]=j2r$warn.p
    nbd.warn.t[i]=j2r$warn.t
    #print(i)
     err[i] <- 0
    }, warning = function(w){return(1)})
  }
 
  # excludes cases with warnings
  co_powerpois<-sum(subset(co_poisson.pval,err==0)<alpha)/length(subset(co_poisson.pval,err==0))
  co_powerqpois<-sum(subset(co_quasipoisson.pval,err==0)<alpha)/length(subset(co_quasipoisson.pval,err==0))
  co_powernb<-sum(subset(co_nbd.pval,err==0)<alpha)/length(subset(co_nbd.pval,err==0))

  dl_powerpois<-sum(subset(dl_poisson.pval,err==0)<alpha)/length(subset(dl_poisson.pval,err==0)) 
  dl_powerqpois<-sum(subset(dl_quasipoisson.pval,err==0)<alpha)/length(subset(dl_quasipoisson.pval,err==0)) 
  dl_powernb<-sum(subset(dl_nbd.pval,err==0)<alpha)/length(subset(dl_nbd.pval,err==0)) 
  
  j2r_powerpois<-sum(subset(j2r_poisson.pval,err==0)<alpha)/length(subset(j2r_poisson.pval,err==0))
  j2r_powerqpois<-sum(subset(j2r_quasipoisson.pval,err==0)<alpha)/length(subset(j2r_quasipoisson.pval,err==0)) 
  j2r_powernb<-sum(subset(j2r_nbd.pval,err==0)<alpha)/length(subset(j2r_nbd.pval,err==0)) 
  
  j2r_powerpoisa<-sum(subset(j2r_poisson.apval,err==0)<alpha)/length(subset(j2r_poisson.apval,err==0))
  j2r_powerqpoisa<-sum(subset(j2r_quasipoisson.apval,err==0)<alpha)/length(subset(j2r_quasipoisson.apval,err==0))
  j2r_powernba<-sum(subset(j2r_nbd.apval,err==0)<alpha)/length(subset(j2r_nbd.apval,err==0))
  
  efy_powerpois<-sum(subset(efy_poisson.pval,err==0)<alpha)/length(subset(efy_poisson.pval,err==0))
  efy_powerqpois<-sum(subset(efy_quasipoisson.pval,err==0)<alpha)/length(subset(efy_quasipoisson.pval,err==0))
  efy_powernb<-sum(subset(efy_nbd.pval,err==0)<alpha)/length(subset(efy_nbd.pval,err==0))
  
  efy_powerpoisa<-sum(subset(efy_poisson.apval,err==0)<alpha)/length(subset(efy_poisson.apval,err==0))
  efy_powerqpoisa<-sum(subset(efy_quasipoisson.apval,err==0)<alpha)/length(subset(efy_quasipoisson.apval,err==0))
  efy_powernba<-sum(subset(efy_nbd.apval,err==0)<alpha)/length(subset(efy_nbd.apval,err==0))
  
  power <- as.data.frame(cbind(co_powerpois,co_powerqpois,co_powernb,
                               dl_powerpois,dl_powerqpois,dl_powernb,
                               j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba))
  
  rratio <- as.data.frame(cbind(co_poisson.rratio,co_quasipoisson.rratio,co_nbd.rratio,
                                dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,
                                j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,
                                efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio))
  
  se <- as.data.frame(cbind(co_poisson.se,co_quasipoisson.se,co_nbd.se,
                            dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,
                            j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,
                            efy_poisson.se,efy_quasipoisson.se,efy_nbd.se,
                            co_poisson.se2,co_quasipoisson.se2,co_nbd.se2,
                            dl_poisson.se2,dl_quasipoisson.se2,dl_nbd.se2,
                            j2r_poisson.se2,j2r_quasipoisson.se2,j2r_nbd.se2,
                            efy_poisson.se2,efy_quasipoisson.se2,efy_nbd.se2))
  
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,
                               efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf))
  
  dl_summnbtheta<-summary(dl_nbd.theta)
  
  dl_nbthetawarn <- sum(dl_nbd.warn)
  nbthetawarn.p <- sum(nbd.warn.p)
  nbthetawarn.t <- sum(nbd.warn.t)
  warnalltype <- sum(err==1)
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
  
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  result <- list(err=err,power=power,dl_theta=dl_summnbtheta,rratio=rratio,se=se,mi_df=mi_df,dl_nbthetawarn=dl_nbthetawarn,dropout=dropout,censortime=censortime,ERgt75p=ERgt75p,ERle75p=ERle75p, py.p=py.p, py.t=py.t, 
                 mr.dl.p=mr.dl.p, mr.ef.p=mr.ef.p, mr.co.p=mr.co.p, mr.j2r.p=mr.j2r.p, events.dl.p=events.dl.p, events.ef.p=events.ef.p, events.j2r.p=events.j2r.p, events.co.p=events.co.p,
                 mr.dl.t=mr.dl.t, mr.ef.t=mr.ef.t, mr.co.t=mr.co.t, mr.j2r.t=mr.j2r.t, events.dl.t=events.dl.t, events.ef.t=events.ef.t, events.j2r.t=events.j2r.t, events.co.t=events.co.t, 
                 p.rate=p.rate, t.rate=t.rate, gamma.p=gamma.p, gamma.t=gamma.t, nbthetawarn.p=nbthetawarn.p, nbthetawarn.t=nbthetawarn.t, warnalltype=warnalltype,outliers.ef=outliers.ef, outliers.j2r=outliers.j2r)
  return(result)
}




