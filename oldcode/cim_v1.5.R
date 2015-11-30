# Missing Data Project Simulation Program
# Author: Robert Wan
# Last Modified: January 21, 2015

# controlled imputation method
cim <- function(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, imp_choice, N.mi, equal_dispersion=1, imp_choice2, delta){  
  
  # Step 1
  if (equal_dispersion==1){
    Y <- c(nb_incomplete.placebo[[1]], nb_incomplete.treatment[[1]])
    X <- c(rep(0,n),rep(1,n))
    censor.time <- c(nb_incomplete.placebo[[2]], nb_incomplete.treatment[[2]])
    nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    gamma.placebo <- nbd.summary$theta
    gamma.treatment <- gamma.placebo
    p1.placebo <- exp(nbd.summary$coeffi[1,1])
    p1.treatment <- exp(nbd.summary$coeffi[1,1]+nbd.summary$coeffi[2,1])
    warn.placebo = (!is.null(nbd.summary$th.warn))
    warn.treatment = (!is.null(nbd.summary$th.warn))
  }
  else{
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
  }
  
  # Step 2
  poisson.trt <- rep(NA,N.mi)
  poisson.se <- rep(NA,N.mi)
  poisson.se2 <- rep(NA,N.mi)
  quasipoisson.trt <- rep(NA,N.mi)
  quasipoisson.se <- rep(NA,N.mi)
  quasipoisson.se2 <- rep(NA,N.mi)
  nbd.trt <- rep(NA,N.mi)
  nbd.mp <- rep(NA,N.mi)
  nbd.mt <- rep(NA,N.mi)
  nbd.th <- rep(NA,N.mi)
  nbd.se <- rep(NA,N.mi)
  nbd.se2 <- rep(NA,N.mi)
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
      T1 <- T-nb_incomplete.placebo[[2]][j]
      p.placebo <- p1.placebo*T1/(gamma.placebo+p1.placebo*T)
      if (T1>0){
        gammaj <- gamma.placebo + nb_incomplete.placebo[[1]][j]
        pj <- p.placebo
        uj <- pj/(1-pj)*gammaj
	      if(imp_choice == 4){
          #delta approach 
          uj = uj*delta[1]
        }
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        imputed <- imputed[[1]] + nb_incomplete.placebo[[2]][j]
        nb_impute.placebo[[j]] <- c(nb_incomplete.placebo[[3]][[j]], imputed)
        if(uj>10) {
          outlier=rbind(outlier, c(round(uj,digits=1), round(gamma.placebo,digits=2),round(T1,digits=3), round(p1.placebo,digits=2), nb_incomplete.placebo[[1]][j], sum(!is.na(imputed[[1]])),imp_choice))
        }
      }
      else{
        nb_impute.placebo[[j]] <- nb_incomplete.placebo[[3]][[j]]
      }
    }
     
    nb_impute.treatment <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-nb_incomplete.treatment[[2]][j]
      p.placebo <- p1.placebo*T1/(gamma.treatment+p1.treatment*nb_incomplete.treatment[[2]][j]+p1.placebo*T1)
      p.treatment <- p1.treatment*T1/(gamma.treatment+p1.treatment*T)
      
      if (T1>0){
        if (imp_choice == 1){
          # J2R
          pj <- p.placebo
        }
        else if(imp_choice == 2){
          # Efficacy estimand
          pj <- p.treatment
        }
        else if (imp_choice == 3){
          # Treatment regimen
          pj <- (p.placebo+p.treatment)/2
        }  
        else if(imp_choice == 4){
          # Used for delta approach
          pj <- p.treatment
        }
        else stop("Imputation choice not implemented yet.")
        gammaj <- gamma.treatment + nb_incomplete.treatment[[1]][j]
        uj <- pj/(1-pj)*gammaj
        if(imp_choice == 4){
          #delta approach 
          uj = uj*delta[2]
        }
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        imputed <- imputed[[1]] + nb_incomplete.treatment[[2]][j]
        nb_impute.treatment[[j]] <- c(nb_incomplete.treatment[[3]][[j]], imputed)
	      if(uj>10) {
          outlier=rbind(outlier, c(round(uj, digits=1), round(gamma.treatment, digits=2), round(T1,digits=3), round(p1.treatment,digits=2), nb_incomplete.treatment[[1]][j], sum(!is.na(imputed[[1]])),imp_choice))
        }
      }
      else{
        nb_impute.treatment[[j]] <- nb_incomplete.treatment[[3]][[j]]
      }
    }
    
   # Step 3
   dropout_time.placebo <- rep(T,n)
   dropout_time.treatment <- rep(T,n)
   new_nb_incomplete.treatment <- nb_incomplete(n,nb_impute.treatment,dropout_time.treatment)
   new_nb_incomplete.placebo   <- nb_incomplete(n,nb_impute.placebo,  dropout_time.placebo)
   
   
   Y <- c(new_nb_incomplete.placebo[[1]], new_nb_incomplete.treatment[[1]])
   X <- c(rep(0,n),rep(1,n))
   censor.time <- c(new_nb_incomplete.placebo[[2]], new_nb_incomplete.treatment[[2]])
   
  if (('all' %in% imp_choice2)|('poisson' %in% imp_choice2)) {  
    poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    poisson.trt[i] <- exp(poisson.summary$coeffi[2,1])
    # Corrected SE formula (original scale)
    poisson.se[i] <- poisson.summary$coeffi[2,2]*exp(poisson.summary$coeffi[2,1])
    # SE of diff in log rates
    poisson.se2[i] <- poisson.summary$coeffi[2,2]
  }
    
  if (('all' %in% imp_choice2)|('quasipoisson' %in% imp_choice2)) {  
    quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
    quasipoisson.trt[i] <- exp(quasipoisson.summary$coeffi[2,1])
    # Corrected SE formula
    quasipoisson.se[i] <- quasipoisson.summary$coeffi[2,2]*exp(quasipoisson.summary$coeffi[2,1])
    quasipoisson.se2[i] <- quasipoisson.summary$coeffi[2,2]
  }
    
  if (('all' %in% imp_choice2)|('nbd' %in% imp_choice2)) {
    nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
    nbd.trt[i] <- exp(nbd.summary$coeffi[2,1])
    
    nbd.mp[i] <- exp(nbd.summary$coeffi[1,1])
    nbd.mt[i] <- exp(nbd.summary$coeffi[1,1] + nbd.summary$coeffi[2,1])
    nbd.th[i] = nbd.summary$theta
    
    # Corrected SE formula
    nbd.se[i] <- nbd.summary$coeffi[2,2]*exp(nbd.summary$coeffi[2,1])  
    nbd.se2[i] <- nbd.summary$coeffi[2,2]
  }  
    count.placebo <- count.placebo + sum(new_nb_incomplete.placebo[[1]])
    count.treatment <- count.treatment + sum(new_nb_incomplete.treatment[[1]])
    
    mr = mr + (1/N.mi)*sum(Y[censor.time>0])/sum(censor.time[censor.time>0])
    ev = ev + sum(Y[censor.time>0])/N.mi
    
    Y.p <- new_nb_incomplete.placebo[[1]]
    censor.time.p <- new_nb_incomplete.placebo[[2]]
    
    Y.t <- new_nb_incomplete.treatment[[1]]
    censor.time.t <- new_nb_incomplete.treatment[[2]]
        
    ev.p = ev.p + sum(Y.p[censor.time.p>0])/N.mi
    ev.t = ev.t + sum(Y.t[censor.time.t>0])/N.mi
  }  
  
    mr.p = mean(nbd.mp)
    mr.t = mean(nbd.mt)
  # Rubin's formula 

  count.placebo <- count.placebo/N.mi
  count.treatment <- count.treatment/N.mi
 
  if (('all' %in% imp_choice2)|('poisson' %in% imp_choice2)) {
    #ATD added: combined estimates for diff in log rates
    mi_poisson.trt2 <- mean(log(poisson.trt))
    mi_poisson.se2 <- sqrt((1+1/N.mi)*var(log(poisson.trt))+mean(poisson.se2^2))
    mi_poisson.df2 <- (N.mi-1)*(1+mean(poisson.se2^2)/((1+1/N.mi)*var(log(poisson.trt))))^2
    # Vhat obs; Note: complete data df=2*n-2
    mi_poisson.vhatobs2 <- (1-(1+1/N.mi)*var(log(poisson.trt))/mi_poisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
    # adjusted degrees of freedom
    mi_poisson.adf2 <- 1/(1/mi_poisson.df2+1/mi_poisson.vhatobs2) 
  }
  
  if (('all' %in% imp_choice2)|('quasipoisson' %in% imp_choice2)) {  
    mi_quasipoisson.trt2 <- mean(log(quasipoisson.trt))
    mi_quasipoisson.se2 <- sqrt((1+1/N.mi)*var(log(quasipoisson.trt))+mean(quasipoisson.se2^2))
    mi_quasipoisson.df2 <- (N.mi-1)*(1+mean(quasipoisson.se2^2)/((1+1/N.mi)*var(log(quasipoisson.trt))))^2
    mi_quasipoisson.vhatobs2 <- (1-(1+1/N.mi)*var(log(quasipoisson.trt))/mi_quasipoisson.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
    mi_quasipoisson.adf2 <- 1/(1/mi_quasipoisson.df2+1/mi_quasipoisson.vhatobs2) 
  }
    
  if (('all' %in% imp_choice2)|('nbd' %in% imp_choice2)) {
    mi_nbd.trt2 <- mean(log(nbd.trt))
    mi_nbd.se2 <- sqrt((1+1/N.mi)*var(log(nbd.trt))+mean(nbd.se2^2))
    mi_nbd.df2 <- (N.mi-1)*(1+mean(nbd.se2^2)/((1+1/N.mi)*var(log(nbd.trt))))^2
    mi_nbd.vhatobs2 <- (1-(1+1/N.mi)*var(log(nbd.trt))/mi_nbd.se2^2)*(2*n-2)*(2*n-2+1)/(2*n-2+3) 
    mi_nbd.adf2 <- 1/(1/mi_nbd.df2+1/mi_nbd.vhatobs2) 
    mi_nbd.th = mean(nbd.th)
  }
 
  ######
  
  #return(c(mi_poisson.trt, mi_poisson.se, mi_poisson.se2=mi_poisson.se2, mi_quasipoisson.trt, mi_quasipoisson.se, mi_quasipoisson.se2, mi_nbd.trt, mi_nbd.se, mi_nbd.se2, count.placebo, count.treatment,
   #        mr=mr, ev=ev, mr.p=mr.p, ev.p=ev.p, mr.t=mr.t, ev.t=ev.t, p.rate=p1.placebo, t.rate=p1.treatment))
  
  #MG Changed to list output according to Annie's code
  out=list(mr=mr,ev=ev,mr.p=mr.p,ev.p=ev.p,mr.t=mr.t,ev.t=ev.t,p.rate=p1.placebo, t.rate=p1.treatment,
          warn.p=warn.placebo, warn.t=warn.treatment, outliers=outlier, gamma.p=gamma.placebo, gamma.t=gamma.treatment, gamma=gamma)
  if (('all' %in% imp_choice2)|('poisson' %in% imp_choice2)) {
    out=append(out, list(mi_poisson.trt2=mi_poisson.trt2,mi_poisson.se2=mi_poisson.se2, mi_poisson.df2=mi_poisson.df2, mi_poisson.adf2=mi_poisson.adf2))
  }
  if (('all' %in% imp_choice2)|('quasipoisson' %in% imp_choice2)) {
    out=append(out, list(mi_quasipoisson.trt2=mi_quasipoisson.trt2,mi_quasipoisson.se2=mi_quasipoisson.se2, mi_quasipoisson.df2=mi_quasipoisson.df2, mi_quasipoisson.adf2=mi_quasipoisson.adf2))
  }
  if (('all' %in% imp_choice2)|('nbd' %in% imp_choice2)) {
    out=append(out, list(mi_nbd.trt2=mi_nbd.trt2,mi_nbd.se2=mi_nbd.se2, mi_nbd.df2=mi_nbd.df2, mi_nbd.adf2=mi_nbd.adf2, mi_nbd.th=mi_nbd.th))
  }
  #return<-list(mi_poisson.trt2=mi_poisson.trt2,mi_poisson.se2=mi_poisson.se2, mi_poisson.df2=mi_poisson.df2, mi_poisson.adf2=mi_poisson.adf2,
  #             mi_quasipoisson.trt2=mi_quasipoisson.trt2,mi_quasipoisson.se2=mi_quasipoisson.se2, mi_quasipoisson.df2=mi_quasipoisson.df2, mi_quasipoisson.adf2=mi_quasipoisson.adf2,
  #             mi_nbd.trt2=mi_nbd.trt2,mi_nbd.se2=mi_nbd.se2, mi_nbd.df2=mi_nbd.df2, mi_nbd.adf2=mi_nbd.adf2,mr=mr,ev=ev,mr.p=mr.p,ev.p=ev.p,mr.t=mr.t,ev.t=ev.t,p.rate=p1.placebo, t.rate=p1.treatment,
  #             warn.p=warn.placebo, warn.t=warn.treatment, outliers=outlier, gamma.p=gamma.placebo, gamma.t=gamma.treatment, gamma=gamma)

}