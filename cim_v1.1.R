# Missing Data Project Simulation Program
# Author: Robert Wan
# Last Modified: August 27, 2014

# controlled imputation method
# Missing Data Project Simulation Program
# Author: Robert Wan
# Last Modified: October 7, 2014

# controlled imputation method
cim <- function(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, imp_choice, N.mi){  
  
  # Step 1
  Y <- nb_incomplete.placebo[[1]]
  censor.time <- nb_incomplete.placebo[[2]]
  
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  
  gamma.placebo <- nbd.summary$theta
  p1.placebo <- exp(nbd.summary$coeffi[1,1])
  warn.placebo = (!is.null(nbd.summary$th.warn))
  
  Y <- nb_incomplete.treatment[[1]]
  censor.time <- nb_incomplete.treatment[[2]]
  
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  
  gamma.treatment <- nbd.summary$theta
  p1.treatment <- exp(nbd.summary$coeffi[1,1])
  warn.treatment = (!is.null(nbd.summary$th.warn))
  
  ###################
  Y.p <- nb_incomplete.placebo[[1]]
  Y.t <- nb_incomplete.treatment[[1]]
  Yimp.p = Y.p
  Yimp.t = Y.t
  Y2 <- c(Y.p, Y.t)
  censor.time2 <- c(nb_incomplete.placebo[[2]],nb_incomplete.treatment[[2]])
  X <- c(rep(0,n),rep(1,n))
  nbd.summary2 <- summary(glm.nb(Y2[censor.time2>0] ~ X[censor.time2>0] + offset(log(censor.time2[censor.time2>0]))))  
  
  m.p=exp(nbd.summary2$coeffi[1,1])
  m.t=exp(nbd.summary2$coeffi[1,1] + nbd.summary2$coeffi[2,1])
  gamma <- nbd.summary2$theta
  
  #gamma.placebo <- gamma
  #p1.placebo <- m.p
  #gamma.treatment <- gamma
  #p1.treatment <- m.t
  #######################
  
  # Step 2
  poisson.trt <- rep(NA,N.mi)
  poisson.se <- rep(NA,N.mi)
  poisson.se2 <- rep(NA,N.mi)
  quasipoisson.trt <- rep(NA,N.mi)
  quasipoisson.se <- rep(NA,N.mi)
  quasipoisson.se2 <- rep(NA,N.mi)
  nbd.trt <- rep(NA,N.mi)
  nbd.se <- rep(NA,N.mi)
  nbd.se2 <- rep(NA,N.mi)
  
  # 
  count.placebo <- 0
  count.treatment <- 0
  mr=0
  ev=0
  mr.p=0
  ev.p=0
  mr.t=0
  ev.t=0
  outlier=c()
  
  for (i in 1:N.mi){
    
    nb_impute.placebo <- array(list(NA),dim=n)
    for (j in 1:n){
       # dropout period
      T1 <- T-nb_incomplete.placebo[[2]][j]
       # See Keene's formula on p.91 of TRISTAN paper. Also step 1 in specs document, p.1.
      p.placebo <- p1.placebo*nb_incomplete.placebo[[2]][j]/(p1.placebo*nb_incomplete.placebo[[2]][j]+gamma.placebo)
      if (T1>0){
        # new theta: old theta (i.e. 1/k) + no. of events prior to dropout (y1) 
        gammaj <- gamma.placebo + nb_incomplete.placebo[[1]][j]
        pj <- p.placebo/(1+p.placebo)
        uj <- pj/(1-pj)*gammaj    
        # time to event from dropout time
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj,digits=1), round(gamma.placebo,digits=2),round(T1,digits=3), round(p1.placebo,digits=2), nb_incomplete.placebo[[1]][j], sum(!is.na(imputed[[1]])),imp_choice))
        }
        # imputed event time during dropout period
        imputed <- imputed[[1]] + nb_incomplete.placebo[[2]][j]
        # new list of event times per id
        nb_impute.placebo[[j]] <- c(nb_incomplete.placebo[[3]][[j]], imputed)
      }
      else{
        nb_impute.placebo[[j]] <- nb_incomplete.placebo[[3]][[j]]
      }
    }
    
    nb_impute.treatment <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-nb_incomplete.treatment[[2]][j]
      p.placebo <- p1.placebo*nb_incomplete.treatment[[2]][j]/(p1.placebo*nb_incomplete.treatment[[2]][j]+gamma.placebo)
      p.treatment <- p1.treatment*nb_incomplete.treatment[[2]][j]/(p1.treatment*nb_incomplete.treatment[[2]][j]+gamma.treatment)
      
      if (T1>0){
        if (imp_choice == 1){
          # J2R
          pj <- p.placebo*(1-p.treatment)/(1-p.placebo*p.treatment)
        }
        else if(imp_choice == 2){
          # Efficacy estimand
          pj <- p.treatment/(1+p.treatment)
        }
        else stop("Imputation choice not implemented yet.")
        gammaj <- gamma.treatment + nb_incomplete.treatment[[1]][j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj, digits=1), round(gamma.treatment, digits=2), round(T1,digits=3), round(p1.treatment,digits=2), nb_incomplete.treatment[[1]][j], sum(!is.na(imputed[[1]])),imp_choice))
        }
        
        imputed <- imputed[[1]] + nb_incomplete.treatment[[2]][j]
        nb_impute.treatment[[j]] <- c(nb_incomplete.treatment[[3]][[j]], imputed)
      }
      else{
        nb_impute.treatment[[j]] <- nb_incomplete.treatment[[3]][[j]]
      }
    }
    # Step 3
    dropout_time.placebo <- rep(T,n)
    dropout_time.treatment <- rep(T,n)
    new_nb_incomplete.placebo <- nb_incomplete(n,nb_impute.placebo,dropout_time.placebo)
    new_nb_incomplete.treatment <- nb_incomplete(n,nb_impute.treatment,dropout_time.treatment)
    
    Y <- c(new_nb_incomplete.placebo[[1]], new_nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time <- c(new_nb_incomplete.placebo[[2]], new_nb_incomplete.treatment[[2]])
    
    poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    poisson.trt[i] <- exp(poisson.summary$coeffi[2,1])
    # Corrected SE formula (original scale)
    poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
    # SE of diff in log rates
    poisson.se2[i] <- poisson.summary$coeffi[2,2]
    
    quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    quasipoisson.trt[i] <- exp(quasipoisson.summary$coeffi[2,1])
    # Corrected SE formula
    quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]
    
    nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    nbd.trt[i] <- exp(nbd.summary$coeffi[2,1])
    # Corrected SE formula
    nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])  
    nbd.se2[i] <- nbd.summary$coeffi[2,2]
    
    count.placebo <- count.placebo + sum(new_nb_incomplete.placebo[[1]])
    count.treatment <- count.treatment + sum(new_nb_incomplete.treatment[[1]])
    
    mr = mr + (1/N.mi)*sum(Y[censor.time>0])/sum(censor.time[censor.time>0])
    ev = ev + sum(Y[censor.time>0])/N.mi
    
    Y.p <- new_nb_incomplete.placebo[[1]]
    censor.time.p <- new_nb_incomplete.placebo[[2]]
    
    Y.t <- new_nb_incomplete.treatment[[1]]
    censor.time.t <- new_nb_incomplete.treatment[[2]]
    
    mr.pi=sum(Y.p[censor.time.p>0])/sum(censor.time.p[censor.time.p>0])
    mr.p = mr.p + (1/N.mi)*mr.pi
    ev.p = ev.p + sum(Y[censor.time.p>0])/N.mi
    
    mr.ti=sum(Y.t[censor.time.t>0])/sum(censor.time.t[censor.time.t>0])
    mr.t = mr.t + (1/N.mi)*mr.ti
    ev.t = ev.t + sum(Y.t[censor.time.t>0])/N.mi
  }  
  
  # Rubin's formula 

  count.placebo <- count.placebo/N.mi
  count.treatment <- count.treatment/N.mi
 
  #ATD added: combined estimates for diff in log rates
  mi_poisson.trt2 <- mean(log(poisson.trt))
  mi_poisson.se2 <- sqrt((1+1/N.mi)*var(log(poisson.trt))+mean(poisson.se2^2)) 
    
  mi_quasipoisson.trt2 <- mean(log(quasipoisson.trt))
  mi_quasipoisson.se2 <- sqrt((1+1/N.mi)*var(log(quasipoisson.trt))+mean(quasipoisson.se2^2))
     
  mi_nbd.trt2 <- mean(log(nbd.trt))
  mi_nbd.se2 <- sqrt((1+1/N.mi)*var(log(nbd.trt))+mean(nbd.se2^2))

  mi_poisson.df2 <- (N.mi-1)*(1+mean(poisson.se2^2)/((1+1/N.mi)*var(log(poisson.trt))))^2
  # Vhat obs; Note: complete data df=2*n-2
  mi_poisson.vhatobs2 <- (1-(1+1/N.mi)*var(log(poisson.trt))/mi_poisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  # adjusted degrees of freedom
  mi_poisson.adf2 <- 1/(1/mi_poisson.df2+1/mi_poisson.vhatobs2) 
  
  mi_quasipoisson.df2 <- (N.mi-1)*(1+mean(quasipoisson.se2^2)/((1+1/N.mi)*var(log(quasipoisson.trt))))^2
  mi_quasipoisson.vhatobs2 <- (1-(1+1/N.mi)*var(log(quasipoisson.trt))/mi_quasipoisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_quasipoisson.adf2 <- 1/(1/mi_quasipoisson.df2+1/mi_quasipoisson.vhatobs2) 
  
  mi_nbd.df2 <- (N.mi-1)*(1+mean(nbd.se2^2)/((1+1/N.mi)*var(log(nbd.trt))))^2
  mi_nbd.vhatobs2 <- (1-(1+1/N.mi)*var(log(nbd.trt))/mi_nbd.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_nbd.adf2 <- 1/(1/mi_nbd.df2+1/mi_nbd.vhatobs2) 

 
  ######
  
  #return(c(mi_poisson.trt, mi_poisson.se, mi_poisson.se2=mi_poisson.se2, mi_quasipoisson.trt, mi_quasipoisson.se, mi_quasipoisson.se2, mi_nbd.trt, mi_nbd.se, mi_nbd.se2, count.placebo, count.treatment,
   #        mr=mr, ev=ev, mr.p=mr.p, ev.p=ev.p, mr.t=mr.t, ev.t=ev.t, p.rate=p1.placebo, t.rate=p1.treatment))
  
  #MG Changed to list output according to Annie's code
  return<-list(mi_poisson.trt2=mi_poisson.trt2,mi_poisson.se2=mi_poisson.se2, mi_poisson.df2=mi_poisson.df2, mi_poisson.adf2=mi_poisson.adf2,
               mi_quasipoisson.trt2=mi_quasipoisson.trt2,mi_quasipoisson.se2=mi_quasipoisson.se2, mi_quasipoisson.df2=mi_quasipoisson.df2, mi_quasipoisson.adf2=mi_quasipoisson.adf2,
               mi_nbd.trt2=mi_nbd.trt2,mi_nbd.se2=mi_nbd.se2, mi_nbd.df2=mi_nbd.df2, mi_nbd.adf2=mi_nbd.adf2,mr=mr,ev=ev,mr.p=mr.p,ev.p=ev.p,mr.t=mr.t,ev.t=ev.t,p.rate=p1.placebo, t.rate=p1.treatment,
               warn.p=warn.placebo, warn.t=warn.treatment, outliers=outlier, gamma.p=gamma.placebo, gamma.t=gamma.treatment, gamma=gamma,p.rate2=m.p, t.rate2=m.t)
}






cim_old <- function(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, imp_choice, N.mi){  
  # Step 1
  Y <- nb_incomplete.placebo[[1]]
  censor.time <- nb_incomplete.placebo[[2]]
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  gamma.placebo <- nbd.summary$theta
  warn.placebo = (!is.null(nbd.summary$th.warn))
# p.placebo <- nbd.summary$theta*exp(nbd.summary$coeffi[1,1])/(1+nbd.summary$theta*exp(nbd.summary$coeffi[1,1]))
# p.placebo <- exp(nbd.summary$coeffi[1,1])/(exp(nbd.summary$coeffi[1,1])+nbd.summary$theta)
  p.placebo <- exp(nbd.summary$coeffi[1,1])*T/(exp(nbd.summary$coeffi[1,1])*T+nbd.summary$theta)
  p.rate <- exp(nbd.summary$coeffi[1,1])
  
  Y <- nb_incomplete.treatment[[1]]
  censor.time <- nb_incomplete.treatment[[2]]
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  gamma.treatment <- nbd.summary$theta
  warn.treatment = (!is.null(nbd.summary$th.warn))
# p.treatment <- nbd.summary$theta*exp(nbd.summary$coeffi[1,1])/(1+nbd.summary$theta*exp(nbd.summary$coeffi[1,1]))
# p.treatment <- exp(nbd.summary$coeffi[1,1])/(exp(nbd.summary$coeffi[1,1])+nbd.summary$theta)
  p.treatment <- exp(nbd.summary$coeffi[1,1])*T/(exp(nbd.summary$coeffi[1,1])*T+nbd.summary$theta)
  t.rate <- exp(nbd.summary$coeffi[1,1])

  # Step 2
  poisson.trt <- rep(NA,N.mi)
  poisson.se <- rep(NA,N.mi)
  poisson.se2 <- rep(NA,N.mi)
  quasipoisson.trt <- rep(NA,N.mi)
  quasipoisson.se <- rep(NA,N.mi)
  quasipoisson.se2 <- rep(NA,N.mi)
  nbd.trt <- rep(NA,N.mi)
  nbd.se <- rep(NA,N.mi)
  nbd.se2 <- rep(NA,N.mi)
  mr=0
  ev=0
  mr.p=0
  ev.p=0
  mr.t=0
  ev.t=0
  outlier=c()
  
  for (i in 1:N.mi){
    nb_impute.placebo <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-nb_incomplete.placebo[[2]][j]
      if (T1>0){
        pj <- p.placebo/(1+p.placebo)
        gammaj <- gamma.placebo + nb_incomplete.placebo[[1]][j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1, n=1, rate=uj, dispersion=1/gammaj)
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj,digits=1), round(gamma.placebo,digits=2),round(T1,digits=3), round(p.rate,digits=2), nb_incomplete.placebo[[1]][j], sum(!is.na(imputed[[1]])),imp_choice))
        }
        
        imputed <- imputed[[1]] + nb_incomplete.placebo[[2]][j]
        nb_impute.placebo[[j]] <- c(nb_incomplete.placebo[[3]][[j]], imputed)
      }
      else{
        nb_impute.placebo[[j]] <- nb_incomplete.placebo[[3]][[j]]
      }
    }
    
    nb_impute.treatment <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-nb_incomplete.treatment[[2]][j]
      if (T1>0){
        if (imp_choice == 1){
          # J2R
          pj <- p.placebo*(1-p.treatment)/(1-p.placebo*p.treatment)
        }
        else if(imp_choice == 2){
          # Efficacy estimand
          pj <- p.treatment/(1+p.treatment)
        }
        else stop("Imputation choice not implemented yet.")
        gammaj <- gamma.treatment + nb_incomplete.treatment[[1]][j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1, n=1, rate=uj, dispersion=1/gammaj)
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj, digits=1), round(gamma.treatment, digits=2), round(T1,digits=3), round(t.rate,digits=2), nb_incomplete.treatment[[1]][j], sum(!is.na(imputed[[1]])),imp_choice))
        }
          
        imputed <- imputed[[1]] + nb_incomplete.treatment[[2]][j]
        nb_impute.treatment[[j]] <- c(nb_incomplete.treatment[[3]][[j]], imputed)
      }
      else{
        nb_impute.treatment[[j]] <- nb_incomplete.treatment[[3]][[j]]
      }
    }
    
   # Step 3
   dropout_time.placebo <- rep(T,n)
   dropout_time.treatment <- rep(T,n)
   new_nb_incomplete.placebo <- nb_incomplete(n,nb_impute.placebo,dropout_time.placebo)
   new_nb_incomplete.treatment <- nb_incomplete(n,nb_impute.treatment,dropout_time.treatment)
   
   Y <- c(new_nb_incomplete.placebo[[1]], new_nb_incomplete.treatment[[1]])
   X <- c(rep(0,n),rep(1,n))
   censor.time <- c(new_nb_incomplete.placebo[[2]], new_nb_incomplete.treatment[[2]])
   
   poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
   poisson.trt[i] <- exp(poisson.summary$coeffi[2,1])
   # Corrected SE formula
   poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
   poisson.se2[i] <- poisson.summary$coeffi[2,2]

   quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
   quasipoisson.trt[i] <- exp(quasipoisson.summary$coeffi[2,1])
   # Corrected SE formula
   quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
   quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]   

   nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
   nbd.trt[i] <- exp(nbd.summary$coeffi[2,1])
   # Corrected SE formula
   nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])
   nbd.se2[i] <- nbd.summary$coeffi[2,2] 
    
    mr = mr + (1/N.mi)*sum(Y[censor.time>0])/sum(censor.time[censor.time>0])
    ev = ev + sum(Y[censor.time>0])/N.mi
    
    Y.p <- new_nb_incomplete.placebo[[1]]
    censor.time.p <- new_nb_incomplete.placebo[[2]]
    
    Y.t <- new_nb_incomplete.treatment[[1]]
    censor.time.t <- new_nb_incomplete.treatment[[2]]
    
    mr.pi=sum(Y.p[censor.time.p>0])/sum(censor.time.p[censor.time.p>0])
    mr.p = mr.p + (1/N.mi)*mr.pi
    ev.p = ev.p + sum(Y[censor.time.p>0])/N.mi
    
    mr.ti=sum(Y.t[censor.time.t>0])/sum(censor.time.t[censor.time.t>0])
    mr.t = mr.t + (1/N.mi)*mr.ti
    ev.t = ev.t + sum(Y.t[censor.time.t>0])/N.mi
    
    #if (mr.pi>3*p.rate) outlier=rbind(outlier, c("Placebo", round(p.rate,digits=3), round(mr.pi,digits=3),round(umax.p,digits=3),round(umax.p,digits=3)))
    #if (mr.ti>3*t.rate) outlier=rbind(outlier, c("Active", round(t.rate,digits=3), round(mr.ti,digits=3),round(umax.t,digits=3)))
    #if (mr.pi>2*p.rate|mr.ti>2*t.rate) print(outlier)
  }  
 
  # Rubin's formula
  mi_poisson.trt <- mean(poisson.trt)
  mi_poisson.se <- sqrt((1+1/N.mi)*var(poisson.trt)+mean(poisson.se^2)) 
  mi_poisson.se2 <- sqrt((1+1/N.mi)*var(log(poisson.trt))+mean(poisson.se2^2))
 
#MG added from Annies code
  mi_poisson.df <- (N.mi-1)*(1+mean(poisson.se2^2)/((1+1/N.mi)*var(poisson.trt)))^2
  # Vhat obs; Note: complete data df=2*n-2
  mi_poisson.vhatobs <- (1-(1+1/N.mi)*var(poisson.trt)/mi_poisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  # adjusted degrees of freedom
  mi_poisson.adf <- 1/(1/mi_poisson.df+1/mi_poisson.vhatobs) 
######
  
  mi_quasipoisson.trt <- mean(quasipoisson.trt)
  mi_quasipoisson.se <- sqrt((1+1/N.mi)*var(quasipoisson.trt)+mean(quasipoisson.se^2))
  mi_quasipoisson.se2 <- sqrt((1+1/N.mi)*var(log(quasipoisson.trt))+mean(quasipoisson.se2^2))


#MG added from Annies code 
  mi_quasipoisson.df <- (N.mi-1)*(1+mean(quasipoisson.se2^2)/((1+1/N.mi)*var(quasipoisson.trt)))^2
  mi_quasipoisson.vhatobs <- (1-(1+1/N.mi)*var(quasipoisson.trt)/mi_quasipoisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_quasipoisson.adf <- 1/(1/mi_quasipoisson.df+1/mi_quasipoisson.vhatobs) 
######

  mi_nbd.trt <- mean(nbd.trt)
  mi_nbd.se <- sqrt((1+1/N.mi)*var(nbd.trt)+mean(nbd.se^2))
  mi_nbd.se2 <- sqrt((1+1/N.mi)*var(log(nbd.trt))+mean(nbd.se2^2))

#MG added from Annies code
  mi_nbd.df <- (N.mi-1)*(1+mean(nbd.se2^2)/((1+1/N.mi)*var(nbd.trt)))^2
  mi_nbd.vhatobs <- (1-(1+1/N.mi)*var(nbd.trt)/mi_nbd.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_nbd.adf <- 1/(1/mi_nbd.df+1/mi_nbd.vhatobs)
######  
  #return(c(mi_poisson.trt, mi_poisson.se, mi_quasipoisson.trt, mi_quasipoisson.se, mi_nbd.trt, mi_nbd.se))

#MG Changed to list output according to Annie's code
 return<-list(mi_poisson.trt=mi_poisson.trt, mi_poisson.se=mi_poisson.se, mi_poisson.se2=mi_poisson.se2, mi_poisson.df=mi_poisson.df, mi_poisson.adf=mi_poisson.adf,
          mi_quasipoisson.trt=mi_quasipoisson.trt, mi_quasipoisson.se=mi_quasipoisson.se, mi_quasipoisson.se2=mi_quasipoisson.se2, mi_quasipoisson.df=mi_quasipoisson.df, mi_quasipoisson.adf=mi_quasipoisson.adf,
          mi_nbd.trt=mi_nbd.trt, mi_nbd.se=mi_nbd.se, mi_nbd.se2=mi_nbd.se2, mi_nbd.df=mi_nbd.df, mi_nbd.adf=mi_nbd.adf,mr=mr,ev=ev,mr.p=mr.p,ev.p=ev.p,mr.t=mr.t,ev.t=ev.t,p.rate=p.rate,t.rate=t.rate, outliers=outlier,
              warn.p=warn.placebo, warn.t=warn.treatment, gamma.p=gamma.placebo, gamma.t=gamma.treatment)
}


# controlled imputation method
cim_imponly <- function(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, imp_choice, N.mi){  
  # Step 1
  Y <- nb_incomplete.placebo[[1]]
  censor.time <- nb_incomplete.placebo[[2]]
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  gamma.placebo <- nbd.summary$theta
  # p.placebo <- nbd.summary$theta*exp(nbd.summary$coeffi[1,1])/(1+nbd.summary$theta*exp(nbd.summary$coeffi[1,1]))
  # p.placebo <- exp(nbd.summary$coeffi[1,1])/(exp(nbd.summary$coeffi[1,1])+nbd.summary$theta)
  p.placebo <- exp(nbd.summary$coeffi[1,1])*T/(exp(nbd.summary$coeffi[1,1])*T+nbd.summary$theta)
  
  Y <- nb_incomplete.treatment[[1]]
  censor.time <- nb_incomplete.treatment[[2]]
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  gamma.treatment <- nbd.summary$theta
  # p.treatment <- nbd.summary$theta*exp(nbd.summary$coeffi[1,1])/(1+nbd.summary$theta*exp(nbd.summary$coeffi[1,1]))
  # p.treatment <- exp(nbd.summary$coeffi[1,1])/(exp(nbd.summary$coeffi[1,1])+nbd.summary$theta)
  p.treatment <- exp(nbd.summary$coeffi[1,1])*T/(exp(nbd.summary$coeffi[1,1])*T+nbd.summary$theta)
  data=data.frame()
  for (i in 1:N.mi){
    nb_impute.placebo <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-nb_incomplete.placebo[[2]][j]
      if (T1>0){
        pj <- p.placebo/(1+p.placebo)
        gammaj <- gamma.placebo + nb_incomplete.placebo[[1]][j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        imputed <- imputed[[1]] + nb_incomplete.placebo[[2]][j]
        nb_impute.placebo[[j]] <- c(nb_incomplete.placebo[[3]][[j]], imputed)
      }
      else{
        nb_impute.placebo[[j]] <- nb_incomplete.placebo[[3]][[j]]
      }
    }
    
    nb_impute.treatment <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-nb_incomplete.treatment[[2]][j]
      if (T1>0){
        if (imp_choice == 1){
          # J2R
          pj <- p.placebo*(1-p.treatment)/(1-p.placebo*p.treatment)
        }
        else if(imp_choice == 2){
          # Efficacy estimand
          pj <- p.treatment/(1+p.treatment)
        }
        else stop("Imputation choice not implemented yet.")
        gammaj <- gamma.treatment + nb_incomplete.treatment[[1]][j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        imputed <- imputed[[1]] + nb_incomplete.treatment[[2]][j]
        nb_impute.treatment[[j]] <- c(nb_incomplete.treatment[[3]][[j]], imputed)
      }
      else{
        nb_impute.treatment[[j]] <- nb_incomplete.treatment[[3]][[j]]
      }
    }
    # Step 3
    dropout_time.placebo <- rep(T,n)
    dropout_time.treatment <- rep(T,n)
    new_nb_incomplete.placebo <- nb_incomplete(n,nb_impute.placebo,dropout_time.placebo)
    new_nb_incomplete.treatment <- nb_incomplete(n,nb_impute.treatment,dropout_time.treatment)
    
    Y.c <- c(nb_complete.placebo[[1]], nb_complete.treatment[[1]])
    Y.i <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    Y.imp <- c(new_nb_incomplete.placebo[[1]], new_nb_incomplete.treatment[[1]])
 
    data=rbind(data, data.frame(imp=i, Imputed=Y.imp))
  }
  return(data)
}  

makeDataFrame = function(nb_incomplete.placebo,nb_incomplete.treatment,nb_complete.placebo,nb_complete.treatment) {
  n=length(nb_incomplete.placebo[[1]])
  Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
  X <- c(rep(0,n),rep(1,n))
  censor.time.p <- mean(subset(nb_incomplete.placebo[[2]],  nb_incomplete.placebo[[2]]<T))
  censor.time.t <- mean(subset(nb_incomplete.treatment[[2]],nb_incomplete.treatment[[2]]<T))
  censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
  
  id=1:n
  Xchar=rep("",length(X)) 
  Xchar[X==0]="Placebo"
  Xchar[X==1]="Active"
  
  data=data.frame(ID=id,TRTn=X,TRT=Xchar,COUNT=Y, STUDY_TIME=censor.time)
  return(data)
}


# controlled imputation method
# Alternative input format - by MG
cimd <- function(T, n, indata, imp_choice, N.mi){  
  
  incmpl.count.p <- indata$COUNT[indata$TRT=="Placebo"]
  incmpl.time.p <- indata$STUDY_TIME[indata$TRT=="Placebo"]  
  incmpl.count.t <- indata$COUNT[indata$TRT=="Active"]
  incmpl.time.t <- indata$STUDY_TIME[indata$TRT=="Active"]

  # Step 1
  Y <- incmpl.count.p
  censor.time <- incmpl.time.p
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  gamma.placebo <- nbd.summary$theta
  warn.placebo = (!is.null(nbd.summary$th.warn))
# p.placebo <- nbd.summary$theta*exp(nbd.summary$coeffi[1,1])/(1+nbd.summary$theta*exp(nbd.summary$coeffi[1,1]))
# p.placebo <- exp(nbd.summary$coeffi[1,1])/(exp(nbd.summary$coeffi[1,1])+nbd.summary$theta)
  p.placebo <- exp(nbd.summary$coeffi[1,1])*T/(exp(nbd.summary$coeffi[1,1])*T+nbd.summary$theta)
  p.rate <- exp(nbd.summary$coeffi[1,1])
  
  Y <- incmpl.count.t
  censor.time <- incmpl.time.t
  nbd.summary <- summary(glm.nb(Y[censor.time>0]~offset(log(censor.time[censor.time>0]))))
  gamma.treatment <- nbd.summary$theta
  warn.treatment = (!is.null(nbd.summary$th.warn))
# p.treatment <- nbd.summary$theta*exp(nbd.summary$coeffi[1,1])/(1+nbd.summary$theta*exp(nbd.summary$coeffi[1,1]))
# p.treatment <- exp(nbd.summary$coeffi[1,1])/(exp(nbd.summary$coeffi[1,1])+nbd.summary$theta)
  p.treatment <- exp(nbd.summary$coeffi[1,1])*T/(exp(nbd.summary$coeffi[1,1])*T+nbd.summary$theta)
  t.rate <- exp(nbd.summary$coeffi[1,1])

  # Step 2
  poisson.trt <- rep(NA,N.mi)
  poisson.se <- rep(NA,N.mi)
  poisson.se2 <- rep(NA,N.mi)
  quasipoisson.trt <- rep(NA,N.mi)
  quasipoisson.se <- rep(NA,N.mi)
  quasipoisson.se2 <- rep(NA,N.mi)
  nbd.trt <- rep(NA,N.mi)
  nbd.se <- rep(NA,N.mi)
  nbd.se2 <- rep(NA,N.mi)
  mr=0
  ev=0
  mr.p=0
  ev.p=0
  mr.t=0
  ev.t=0
  outlier=c()
  for (i in 1:N.mi){
    umax.p=0
    nb_impute.placebo <- rep(NA,n)
    for (j in 1:n){
      T1 <- T-incmpl.time.p[j]
      if (T1>0){
        pj <- p.placebo/(1+p.placebo)
        gammaj <- gamma.placebo + incmpl.count.p[j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1, n=1, rate=uj/T1, dispersion=1/gammaj)
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj,digits=1), round(gamma.placebo,digits=2),round(T1,digits=3), round(p.rate,digits=2), incmpl.count.p[j], sum(!is.na(imputed[[1]])),imp_choice))
        }
        
	  nb_impute.placebo[j] <- incmpl.count.p[j] + sum(!is.na(imputed[[1]]))
      }
      else{
        nb_impute.placebo[j] <- incmpl.count.p[j]
      }
    }
    
    nb_impute.treatment <- rep(NA,n)
    umax.t=0
    for (j in 1:n){
      T1 <- T-incmpl.time.p[j]
      if (T1>0){
        if (imp_choice == 1){
          # J2R
          pj <- p.placebo*(1-p.treatment)/(1-p.placebo*p.treatment)
        }
        else if(imp_choice == 2){
          # Efficacy estimand
          pj <- p.treatment/(1+p.treatment)
        }
        else stop("Imputation choice not implemented yet.")
        gammaj <- gamma.treatment + incmpl.count.p[j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1, n=1, rate=uj/T1, dispersion=1/gammaj)
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj, digits=1), round(gamma.treatment, digits=2), round(T1,digits=3), round(t.rate,digits=2), incmpl.count.t[j], sum(!is.na(imputed[[1]])),imp_choice))
        }
          
        nb_impute.treatment[j] <- incmpl.count.t[j] + sum(!is.na(imputed[[1]]))
	print(c(incmpl.count.t[j], sum(!is.na(imputed[[1]]))))
      }
      else{
        nb_impute.treatment[j] <- incmpl.count.t[j]
      }
    }
    
   # Step 3
   dropout_time.placebo <- rep(T,n)
   dropout_time.treatment <- rep(T,n)
 
   Y <- c(nb_impute.placebo, nb_impute.treatment)
   X <- c(rep(0,n),rep(1,n))
   censor.time <- c(dropout_time.placebo, dropout_time.treatment)
   
   poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
   poisson.trt[i] <- exp(poisson.summary$coeffi[2,1])
   # Corrected SE formula
   poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
   poisson.se2[i] <- poisson.summary$coeffi[2,2]

   quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
   quasipoisson.trt[i] <- exp(quasipoisson.summary$coeffi[2,1])
   # Corrected SE formula
   quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
   quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]   

   nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
   nbd.trt[i] <- exp(nbd.summary$coeffi[2,1])
   # Corrected SE formula
   nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])
   nbd.se2[i] <- nbd.summary$coeffi[2,2] 
    
    mr = mr + (1/N.mi)*sum(Y[censor.time>0])/sum(censor.time[censor.time>0])
    ev = ev + sum(Y[censor.time>0])/N.mi
    
    Y.p <- nb_impute.placebo
    censor.time.p <- dropout_time.placebo
    
    Y.t <- nb_impute.treatment
    censor.time.t <- dropout_time.treatment
    
    mr.pi=sum(Y.p[censor.time.p>0])/sum(censor.time.p[censor.time.p>0])
    mr.p = mr.p + (1/N.mi)*mr.pi
    ev.p = ev.p + sum(Y[censor.time.p>0])/N.mi
    
    mr.ti=sum(Y.t[censor.time.t>0])/sum(censor.time.t[censor.time.t>0])
    mr.t = mr.t + (1/N.mi)*mr.ti
    ev.t = ev.t + sum(Y.t[censor.time.t>0])/N.mi
    
    #if (mr.pi>3*p.rate) outlier=rbind(outlier, c("Placebo", round(p.rate,digits=3), round(mr.pi,digits=3),round(umax.p,digits=3),round(umax.p,digits=3)))
    #if (mr.ti>3*t.rate) outlier=rbind(outlier, c("Active", round(t.rate,digits=3), round(mr.ti,digits=3),round(umax.t,digits=3)))
    #if (mr.pi>2*p.rate|mr.ti>2*t.rate) print(outlier)
  }  

  # Rubin's formula
  mi_poisson.trt <- mean(poisson.trt)
  mi_poisson.se <- sqrt((1+1/N.mi)*var(poisson.trt) + mean(poisson.se^2)) 
  mi_poisson.se2 <- sqrt((1+1/N.mi)*var(log(poisson.trt)) + mean(poisson.se2^2))
 
#MG added from Annies code
  mi_poisson.df <- (N.mi-1)*(1+mean(poisson.se2^2)/((1+1/N.mi)*var(poisson.trt)))^2
  # Vhat obs; Note: complete data df=2*n-2
  mi_poisson.vhatobs <- (1-(1+1/N.mi)*var(poisson.trt)/mi_poisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  # adjusted degrees of freedom
  mi_poisson.adf <- 1/(1/mi_poisson.df+1/mi_poisson.vhatobs) 
######
  
  mi_quasipoisson.trt <- mean(quasipoisson.trt)
  mi_quasipoisson.se <- sqrt((1+1/N.mi)*var(quasipoisson.trt)+mean(quasipoisson.se^2))
  mi_quasipoisson.se2 <- sqrt((1+1/N.mi)*var(log(quasipoisson.trt))+mean(quasipoisson.se2^2))


#MG added from Annies code 
  mi_quasipoisson.df <- (N.mi-1)*(1+mean(quasipoisson.se2^2)/((1+1/N.mi)*var(quasipoisson.trt)))^2
  mi_quasipoisson.vhatobs <- (1-(1+1/N.mi)*var(quasipoisson.trt)/mi_quasipoisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_quasipoisson.adf <- 1/(1/mi_quasipoisson.df+1/mi_quasipoisson.vhatobs) 
######

  mi_nbd.trt <- mean(nbd.trt)
  mi_nbd.se <- sqrt((1+1/N.mi)*var(nbd.trt)+mean(nbd.se^2))
  mi_nbd.se2 <- sqrt((1+1/N.mi)*var(log(nbd.trt))+mean(nbd.se2^2))

#MG added from Annies code
  mi_nbd.df <- (N.mi-1)*(1+mean(nbd.se2^2)/((1+1/N.mi)*var(nbd.trt)))^2
  mi_nbd.vhatobs <- (1-(1+1/N.mi)*var(nbd.trt)/mi_nbd.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_nbd.adf <- 1/(1/mi_nbd.df+1/mi_nbd.vhatobs)
######  
  #return(c(mi_poisson.trt, mi_poisson.se, mi_quasipoisson.trt, mi_quasipoisson.se, mi_nbd.trt, mi_nbd.se))

#MG Changed to list output according to Annie's code
 return<-list(mi_poisson.trt=mi_poisson.trt, mi_poisson.se=mi_poisson.se, mi_poisson.se2=mi_poisson.se2, mi_poisson.df=mi_poisson.df, mi_poisson.adf=mi_poisson.adf,
          mi_quasipoisson.trt=mi_quasipoisson.trt, mi_quasipoisson.se=mi_quasipoisson.se, mi_quasipoisson.se2=mi_quasipoisson.se2, mi_quasipoisson.df=mi_quasipoisson.df, mi_quasipoisson.adf=mi_quasipoisson.adf,
          mi_nbd.trt=mi_nbd.trt, mi_nbd.se=mi_nbd.se, mi_nbd.se2=mi_nbd.se2, mi_nbd.df=mi_nbd.df, mi_nbd.adf=mi_nbd.adf,mr=mr,ev=ev,mr.p=mr.p,ev.p=ev.p,mr.t=mr.t,ev.t=ev.t,p.rate=p.rate,t.rate=t.rate, outliers=outlier,
              warn.p=warn.placebo, warn.t=warn.treatment)
}


    

# controlled imputation method
# Following James Rogers approach - by MG
cim_mg <- function(Tm, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, imp_choice, N.mi){  
  # Step 1
  Y.p <- nb_incomplete.placebo[[1]]
  Y.t <- nb_incomplete.treatment[[1]]
  Yimp.p = Y.p
  Yimp.t = Y.t
  Y2 <- c(Y.p, Y.t)
  censor.time2 <- c(nb_incomplete.placebo[[2]],nb_incomplete.treatment[[2]])
  X <- c(rep(0,n),rep(1,n))
  nbd.summary2 <- summary(glm.nb(Y2[censor.time2>0] ~ X[censor.time2>0] + offset(log(censor.time2[censor.time2>0]))))  
  
  m.p=exp(nbd.summary2$coeffi[1,1])
  m.t=exp(nbd.summary2$coeffi[1,1] + nbd.summary2$coeffi[2,1])
  gamma <- nbd.summary2$theta
  warn = (!is.null(nbd.summary2$th.warn))
  
  # Step 2
  poisson.trt <- rep(NA,N.mi)
  poisson.se <- rep(NA,N.mi)
  poisson.se2 <- rep(NA,N.mi)
  quasipoisson.trt <- rep(NA,N.mi)
  quasipoisson.se <- rep(NA,N.mi)
  quasipoisson.se2 <- rep(NA,N.mi)
  nbd.trt <- rep(NA,N.mi)
  nbd.se <- rep(NA,N.mi)
  nbd.se2 <- rep(NA,N.mi)
  mr=0
  ev=0
  mr.p=0
  ev.p=0
  mr.t=0
  ev.t=0
  outlier=c()
  for (i in 1:N.mi){
    for (j in 1:n){
      Tdo=nb_incomplete.treatment[[2]][j]
      T1 <- Tm-Tdo
      if (T1>0){
        gammaj <- gamma + Y.p[j]
        m1=m.p*Tdo
        m2=m.p*T1
        pstar=m2/(gamma+m1+m2)
        pj=pstar/(1-pstar)	
        uj <- pj*gammaj
        imputed <- rnbinom(1,size=gammaj, prob=1-pstar)
        Yimp.p[j]=Y.p[j]+imputed
        
        print(c(Tdo,T1,m1,m2,pstar,Yimp.p[j],Y.p[j],imputed))
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(1-pstar,digits=1), round(gamma,digits=3),round(T1,digits=3), round(m.p,digits=2), Y.p[j], Yimp.p[j],imp_choice))
        }	
      }
    }
    
    for (j in 1:n){
      Tdo=nb_incomplete.treatment[[2]][j]
      T1 <- Tm-Tdo
      if (T1>0){
        if (imp_choice == 1){
          # J2R
          m1=m.t*Tdo
          m2=m.p*T1
        }
        else if(imp_choice == 2){
          # Efficacy estimand
          m1=m.t*Tdo
          m2=m.t*T1	
        }
        else if(imp_choice == 3){
          # CR
          m1=m.p*Tdo
          m2=m.p*T1  
        }
        else stop("Imputation choice not implemented yet.")
        gammaj <- gamma + Y.p[j]
        pstar=m2/(gamma+m1+m2)
        pj=pstar/(1-pstar)	
        uj <- pj*gammaj
        imputed <- rnbinom(1,size=gammaj, prob=1-pstar)
        Yimp.t[j]=Y.t[j]+imputed
        
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj,digits=1), round(gamma,digits=3),round(T1,digits=3), round(m.t,digits=2), Y.t[j], Yimp.t[j],imp_choice))
        }
      }
    }
    
    # Step 3
    
    Y <- c(Yimp.p, Yimp.t)
    X <- c(rep(0,n),rep(1,n))
    censor.time <- c(rep(Tm,n),rep(Tm,n))
    
    poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    poisson.trt[i] <- exp(poisson.summary$coeffi[2,1])
    # Corrected SE formula
    poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
    poisson.se2[i] <- poisson.summary$coeffi[2,2]
    
    quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    quasipoisson.trt[i] <- exp(quasipoisson.summary$coeffi[2,1])
    # Corrected SE formula
    quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]   
    
    nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    nbd.trt[i] <- exp(nbd.summary$coeffi[2,1])
    # Corrected SE formula
    nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])
    nbd.se2[i] <- nbd.summary$coeffi[2,2] 
    
    mr = mr + (1/N.mi)*sum(Y[censor.time>0])/sum(censor.time[censor.time>0])
    ev = ev + sum(Y[censor.time>0])/N.mi
    
    mr.pi=sum(Yimp.p)/(Tm*n)
    mr.p = mr.p + (1/N.mi)*mr.pi
    ev.p = ev.p + sum(Y.p)/N.mi
    
    mr.ti=sum(Yimp.t)/(Tm*n)
    mr.t = mr.t + (1/N.mi)*mr.ti
    ev.t = ev.t + sum(Y.t)/N.mi
  }  
  
  #print(c(m.p,mr.p,m.t,mr.t), digits=3)
  
  # Rubin's formula
  mi_poisson.trt <- mean(poisson.trt)
  mi_poisson.se <- sqrt((1+1/N.mi)*var(poisson.trt)+mean(poisson.se^2)) 
  mi_poisson.se2 <- sqrt((1+1/N.mi)*var(log(poisson.trt))+mean(poisson.se2^2))
  
  #MG added from Annies code
  mi_poisson.df <- (N.mi-1)*(1+mean(poisson.se2^2)/((1+1/N.mi)*var(poisson.trt)))^2
  # Vhat obs; Note: complete data df=2*n-2
  mi_poisson.vhatobs <- (1-(1+1/N.mi)*var(poisson.trt)/mi_poisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  # adjusted degrees of freedom
  mi_poisson.adf <- 1/(1/mi_poisson.df+1/mi_poisson.vhatobs) 
  ######
  
  mi_quasipoisson.trt <- mean(quasipoisson.trt)
  mi_quasipoisson.se <- sqrt((1+1/N.mi)*var(quasipoisson.trt)+mean(quasipoisson.se^2))
  mi_quasipoisson.se2 <- sqrt((1+1/N.mi)*var(log(quasipoisson.trt))+mean(quasipoisson.se2^2))
  
  
  #MG added from Annies code 
  mi_quasipoisson.df <- (N.mi-1)*(1+mean(quasipoisson.se2^2)/((1+1/N.mi)*var(quasipoisson.trt)))^2
  mi_quasipoisson.vhatobs <- (1-(1+1/N.mi)*var(quasipoisson.trt)/mi_quasipoisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_quasipoisson.adf <- 1/(1/mi_quasipoisson.df+1/mi_quasipoisson.vhatobs) 
  ######
  
  mi_nbd.trt <- mean(nbd.trt)
  mi_nbd.se <- sqrt((1+1/N.mi)*var(nbd.trt)+mean(nbd.se^2))
  mi_nbd.se2 <- sqrt((1+1/N.mi)*var(log(nbd.trt))+mean(nbd.se2^2))
  
  #MG added from Annies code
  mi_nbd.df <- (N.mi-1)*(1+mean(nbd.se2^2)/((1+1/N.mi)*var(nbd.trt)))^2
  mi_nbd.vhatobs <- (1-(1+1/N.mi)*var(nbd.trt)/mi_nbd.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
  mi_nbd.adf <- 1/(1/mi_nbd.df+1/mi_nbd.vhatobs)
  ######  
  #return(c(mi_poisson.trt, mi_poisson.se, mi_quasipoisson.trt, mi_quasipoisson.se, mi_nbd.trt, mi_nbd.se))
  
  #MG Changed to list output according to Annie's code
  return<-list(mi_poisson.trt=mi_poisson.trt, mi_poisson.se=mi_poisson.se, mi_poisson.se2=mi_poisson.se2, mi_poisson.df=mi_poisson.df, mi_poisson.adf=mi_poisson.adf,
               mi_quasipoisson.trt=mi_quasipoisson.trt, mi_quasipoisson.se=mi_quasipoisson.se, mi_quasipoisson.se2=mi_quasipoisson.se2, mi_quasipoisson.df=mi_quasipoisson.df, mi_quasipoisson.adf=mi_quasipoisson.adf,
               mi_nbd.trt=mi_nbd.trt, mi_nbd.se=mi_nbd.se, mi_nbd.se2=mi_nbd.se2, mi_nbd.df=mi_nbd.df, mi_nbd.adf=mi_nbd.adf,mr=mr,ev=ev,mr.p=mr.p,ev.p=ev.p,mr.t=mr.t,ev.t=ev.t,p.rate=m.p,t.rate=m.t, outliers=outlier,
               warn.p=warn, warn.t=warn)
}





