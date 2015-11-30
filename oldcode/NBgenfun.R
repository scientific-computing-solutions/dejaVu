
# Complete data generation, nb_complete contains complete event time for each subject

nb_complete <- function(T,n,rate,dispersion){
  t <- array(list(NA),dim=n)
  
  mean <- rate*T
  variance <- mean + dispersion*mean^2
  p <- 1 - mean/variance
  # alt: p <- (dispersion*mean)/(1+dispersion*mean)
  r <- mean*(1-p)/p
  # alt: r <- 1/dispersion
  
  if(dispersion==0){
    poisson.rate <- rep(rate,n)
  } else {
    # note: (1-p)/p is rate, i.e. scale = p/(1-p)
    lambda <- rgamma(n,r,(1-p)/p)
    poisson.rate <- lambda/T
    # alt: lambda <- rgamma(n,r,1)
    # poisson.rate <- mean*(lambda)*dispersion
  }
  
  for (i in 1:n){
    # does not take into account the requirement for the adjacent exacerbations to be > 7 days apart
    wait.time <- rexp(1,poisson.rate[i])
    while (sum(wait.time)<=T){
      wait.time <- c(wait.time,rexp(1,poisson.rate[i]))
    }
    
    if (length(wait.time)==1) {wait.time<-NA}
    else {wait.time <- wait.time[1:(length(wait.time)-1)]}
    
    t[[i]][1] <- wait.time[1]
    if (length(wait.time)>1){
      for (j in seq(2,length(wait.time),by=1)){
        t[[i]] <- append(t[[i]], wait.time[j]+t[[i]][j-1])
      }
    }
  }
  return(t)
}


# Incomplete data from the complete data with dropout information
nb_incomplete <- function(n,nb_complete,dropout_time){
  count <- rep(0,n)
  t <- array(list(NA),dim=n)
  
  for (i in 1:n){
    j <-1
    while ((j<=length(nb_complete[[i]])) && (!is.na(nb_complete[[i]][j])) && (nb_complete[[i]][j]<=dropout_time[i])){
      count[i] <- count[i] +1
      if (j==1){
        t[[i]][1] <- nb_complete[[i]][j]
      }
      else{
        t[[i]] <- append(t[[i]], nb_complete[[i]][j])
      }
      j <-j+1
    }
  }
  return(list(count,dropout_time,t))
}


# 4. MAR with rate variations
dropout_time4 <- function(T,n,events_time,rates,var){
  t <- rep(NA,n)
  for (i in 1:n){
    j <- 1
    while ((j<= length(events_time[[i]])+1) && is.na(t[i])){
      t1 <- rexp(1,rates[j]*exp(rnorm(1,mean=0,sd=sqrt(var))))
      if (j==1){
        end <- 0
      }
      else{
        end <- events_time[[i]][j-1]
      }
      if ((is.na(events_time[[i]][j])) || (j==1+length(events_time[[i]]))){
        t[i] <- pmin(t1+end,T)
      }
      else{
        if (t1+end<=events_time[[i]][j]) {
          t[i] <- t1+end
        }
      }
      j <- j+1
    }
    if (is.na(t[i])) t[i] <- T
  }
  return(t)
}


# Incomplete data based on CR

nb_incompleteCR <- function(marlst,T,rates,var){
    mard=marlst$dat
    mard$COUNT_CR=mard$COUNT
    mard$ST_CR=mard$STUDY_TIME
    mard$cCR=mard$cc
    n=length(mard$TRT)/2
    nbd.summary <- summary(glm.nb(data=mard, COUNT~TRTn + offset(log(STUDY_TIME))))
    gamma <- nbd.summary$theta
    p1.placebo <- exp(nbd.summary$coeffi[1,1])
    p1.treatment <- exp(nbd.summary$coeffi[1,1]+nbd.summary$coeffi[2,1])
    nb_impute.placebo <- array(list(NA),dim=n)
    #Impute using CR method
    nb_impute.treatment <- array(list(NA),dim=n)
    for (j in 1:n){
      T1 <- T-mard$STUDY_TIME[n+j]
      p.placebo <- p1.placebo*T1/(gamma+p1.treatment*mard$STUDY_TIME[n+j]+p1.placebo*T1)
      p.treatment <- p1.treatment*T1/(gamma+p1.treatment*T)
      if (T1>0){
        # Treatment regimen
        pj <- (p.placebo+p.treatment)/2
        gammaj <- gamma + mard$COUNT[j]
        uj <- pj/(1-pj)*gammaj
        imputed <- nb_complete(T=T1,n=1,rate=uj/T1,dispersion=1/gammaj)
        imputed <- imputed[[1]] + mard$STUDY_TIME[n+j]
        nb_impute.treatment[[j]] <- c(marlst$times.t[[j]], imputed)
      }
      else{
        nb_impute.treatment[[j]] <- marlst$times.t[[j]]
      }
    }
    
    newdo = dropout_time4(T,n,nb_impute.treatment,rates,var)
    newincompl = nb_incomplete(n,nb_impute.treatment, newdo)
    mard$COUNT_CR[(n+1):(2*n)]=newincompl[[1]]
    cc.t=rep(0,n)
    for (k in 1:n) cc.t[k]=sum(!is.na(nb_impute.treatment[[k]]))
    mard$cCR[(n+1):(2*n)]=cc.t
    mard$ST_CR[(n+1):(2*n)]=newdo
    
    return(list(dat=mard,times.p=marlst$times.p, times.t=marlst$times.t, times.t.cr=newincompl[[3]]))
}

#tstd=genNmar(N=1,n=1000,rate.placebo=0.9,rr=0.1,T=1,
#             dispersion.placebo=1.1,dispersion.treatment=1.1,
#             rates=seq(0.001,20,1),
#             drop_rate_variance=0.05)
#tstdcr=nb_incompleteCR(tstd,1,seq(0.001,20,1),0.05)
#summary(tstdcr$dat)

