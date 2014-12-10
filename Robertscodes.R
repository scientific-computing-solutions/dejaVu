
# Complete data generation, nb_complete contains complete event time for each subject
##' @T 
##' @n 
##' @rate 
##' @dispersion 
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

# Incomplete data from the complete data with dropout information - old
#nb_incomplete <- function(n,nb_complete,dropout_time){
#  count <- rep(0,n)
#  for (i in 1:n){
#    j <-1
#    while ((j<=length(nb_complete[[i]])) && (!is.na(nb_complete[[i]][j])) && (nb_complete[[i]][j]<=dropout_time[i])){
#      count[i] <- count[i] +1
#      j <-j+1
#    }
#  }
#  return(list(count,dropout_time))
#}


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


# 1. MCAR with constant rate
dropout_time1 <- function(T,n,rate){
  t <- pmin(rexp(n,rate),T)
  return(t)
}

# 2. MCAR with rate variations
dropout_time2 <- function(T,n,rate,var){
  t <- rep(NA,n)
  for (i in 1:n){
    rate1 <- rate*exp(rnorm(1,mean=0,sd=sqrt(var)))
    t[i] <- rexp(1,rate1)
  }
  t <- pmin(t,T)
  return(t)
}


# 3. MAR with constant rates, rates is a vector whose length is determined by the maximum number of events
dropout_time3 <- function(T,n,events_time,rates){
  t <- rep(NA,n)
  for (i in 1:n){
    j <- 1
    while ((j<= length(events_time[[i]])+1) && is.na(t[i])){
      # dropout rate increases with observed number of exacerbations, i.e. t1 decreases with increasing # of exacerbations  
      t1 <- rexp(1,rates[j])
      if (j==1){
        end <- 0
      }
      else{
        end <- events_time[[i]][j-1]
      }
      if ((is.na(events_time[[i]][j])) | (j==1+length(events_time[[i]]))){
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

