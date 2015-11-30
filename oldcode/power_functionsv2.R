# Missing Data Project Simulation Program
# 

library(MASS)
library(plyr)

## libraries to paralellise the loops
#library(foreach)
#library(doMC)



makeDataFrame = function(nb_incomplete.placebo,nb_incomplete.treatment) {
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


genNmar=function(N,n,rate.placebo,rr,T,dispersion.placebo,dispersion.treatment,rates,drop_rate_variance){
  data=data.frame()
  tim.p=c()
  tim.t=c()
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
    
    cc.p=rep(0,n)
    for (k in 1:n) cc.p[k]=sum(!is.na(nb_complete.placebo[[k]]))
    cc.t=rep(0,n)
    for (k in 1:n) cc.t[k]=sum(!is.na(nb_complete.treatment[[k]]))
    
    id=1:n
    Xchar=rep("",length(X)) 
    Xchar[X==0]="Placebo"
    Xchar[X==1]="Treatment"
    
    data=rbind(data, data.frame(sim=i, ID=id, TRTn=X, TRT=Xchar, COUNT=Y, STUDY_TIME=censor.time, cc=c(cc.p,cc.t) ))
    tim.p=c(tim.p,nb_incomplete.placebo[[3]])
    tim.t=c(tim.t,nb_incomplete.treatment[[3]])
  }
  return(list(dat=data,times.p=tim.p,times.t=tim.t))
}


setdf = function(N.rep) {
    emptyframe=data.frame(rratio=rep(NA,N.rep),se=rep(NA,N.rep),se2=rep(NA,N.rep),
                        theta=rep(NA,N.rep),pval=rep(NA,N.rep),warn=rep(NA,N.rep),
                        apval=rep(NA,N.rep),df=rep(NA,N.rep),adf=rep(NA,N.rep))
    emptyframe2=data.frame(dropout.p = rep(NA,N.rep), dropout.t = rep(NA,N.rep), err= rep(NA,N.rep),
                           events.co.p = rep(NA,N.rep), events.co.t = rep(NA,N.rep),
                           events.dl.p = rep(NA,N.rep), events.dl.t = rep(NA,N.rep),
                           events.ef.p = rep(NA,N.rep), events.ef.t = rep(NA,N.rep),
                           events.j2r.p = rep(NA,N.rep), events.j2r.t = rep(NA,N.rep),
                           events.delta.p = rep(NA,N.rep), events.delta.t = rep(NA,N.rep),
                           events.tr.p = rep(NA,N.rep), events.tr.t = rep(NA,N.rep),
                           mr.co.p = rep(NA,N.rep), mr.co.t = rep(NA,N.rep),
                           mr.dl.p = rep(NA,N.rep), mr.dl.t = rep(NA,N.rep),
                           mr.ef.p = rep(NA,N.rep), mr.ef.t = rep(NA,N.rep),
                           mr.j2r.p = rep(NA,N.rep), mr.j2r.t = rep(NA,N.rep),
                           mr.delta.p = rep(NA,N.rep), mr.delta.t = rep(NA,N.rep),
                           mr.tr.p = rep(NA,N.rep), mr.tr.t = rep(NA,N.rep),
                           j2r_nbd.th  = rep(NA,N.rep), delta_nbd.th  = rep(NA,N.rep),
                           censor.time.p = rep(NA,N.rep), censor.time.t = rep(NA,N.rep),
                           censor.timeT.p = rep(NA,N.rep), censor.timeT.t = rep(NA,N.rep),
                           t.rate = rep(NA,N.rep), p.rate = rep(NA,N.rep),
                           ERgt75p.p = rep(NA,N.rep), ERgt75p.t = rep(NA,N.rep), ERle75p.p = rep(NA,N.rep), ERle75p.t = rep(NA,N.rep))
  rl=data.frame(
  #Complete data
    co_poisson=emptyframe, co_quasipoisson=emptyframe, co_nbd=emptyframe,
  # no imputation (direct likelihood approach)
    dl_poisson=emptyframe, dl_quasipoisson=emptyframe, dl_nbd=emptyframe,
  # jump to reference 
    j2r_poisson=emptyframe, j2r_quasipoisson=emptyframe, j2r_nbd=emptyframe,
  # efficacy estimand
    efy_poisson=emptyframe, efy_quasipoisson=emptyframe, efy_nbd=emptyframe,
  # Treatment regimen estimand
    tr_poisson=emptyframe, tr_quasipoisson=emptyframe, tr_nbd=emptyframe,
  # delta method
    delta_poisson=emptyframe, delta_quasipoisson=emptyframe, delta_nbd=emptyframe,  
  # rates, etc  
    emptyframe2
  )
  return(rl)
}

# MAR with variation in drop-out rate
# 
gc=function(cdat, j, T, n, N.mi, onepertrt, X, rl, meth, meth2, delta ) {
  #rate.placebo, dispersion.placebo, rate.treatment, dispersion.treatment, rates, drop_rate_variance
    rl$err[i] <-tryCatch({
    
    i=1    
    cdati=cdat$dat[cdat$dat$sim==j,] 
    cdi.p=cdati[cdati$TRTn==0,]
    cdi.t=cdati[cdati$TRTn==1,]
    
    # sim ID TRTn     TRT COUNT STUDY_TIME    
       
    rl$dropout.p[i] <- length(cdi.p$STUDY_TIME[cdi.p$STUDY_TIME < T])/n  
    rl$dropout.t[i] <- length(cdi.t$STUDY_TIME[cdi.t$STUDY_TIME < T])/n
    
    nb_inc.p <- cdi.p$COUNT
    nb_inc.t <- cdi.t$COUNT
    do.p <- cdi.p$STUDY_TIME
    do.t <- cdi.t$STUDY_TIME
    nb_complete.placebo=cdi.p$cc
    nb_complete.treatment=cdi.t$cc
      
    # mean ER of patients who stayed for at least 75% of total follow-up period vs early dropouts
    rl$ERgt75p.p[i] <- mean(subset(nb_inc.p/do.p, do.p>0.75*T))
    rl$ERgt75p.t[i] <- mean(subset(nb_inc.t/do.t, do.t>0.75*T))
    rl$ERle75p.p[i] <- mean(subset(nb_inc.p/do.p, do.p<=0.75*T))
    rl$ERle75p.t[i] <- mean(subset(nb_inc.t/do.t, do.t<=0.75*T))
   
    Y <- cdati$COUNT
    censor.time <- cdati$STUDY_TIME
    nb_incomplete.placebo=list(nb_inc.p, do.p, cdat$times.p[(j-1)*n + 1:n])
    nb_incomplete.treatment=list(nb_inc.t, do.t, cdat$times.t[(j-1)*n + 1:n])
    
    rl$censor.time.p[i] <- mean(subset(do.p, do.p<T))
    rl$censor.time.t[i] <- mean(subset(do.t, do.t<T))
    rl$censor.timeT.p[i] <- mean(do.p)
    rl$censor.timeT.t[i] <- mean(do.t)
          
    rl$events.co.p[i]=sum(cdi.p$cc)            
    rl$events.co.t[i]=sum(cdi.t$cc)
    
    rl$events.dl.p[i]=sum(nb_inc.p[do.p>0])
    rl$events.dl.t[i]=sum(nb_inc.t[do.t>0])
    
    #Complete data    
    Yco <- c(cdi.p$cc, cdi.t$cc)
    if (('all' %in% meth2)|('complete' %in% meth2)) {
      if (('all' %in% meth)|('poisson' %in% meth)) {
        co_poisson.summary <- summary(glm(Yco~X+offset(log(rep(T,2*n))), family=poisson()))
      
        #Roberts formula for s.e.
        rl$co_poisson.se[i] <- co_poisson.summary$coeffi[2,2]*exp(co_poisson.summary$coeffi[2,1])
        rl$co_poisson.se2[i] <- co_poisson.summary$coeffi[2,2]
        rl$co_poisson.rratio[i] <- exp(co_poisson.summary$coeffi[2,1])   #Rate ratio
        rl$co_poisson.pval[i] <- co_poisson.summary$coef[2,4] #p-value
      }
    
      if (('all' %in% meth)|('quasipoisson' %in% meth)) {
        co_quasipoisson.summary <- summary(glm(Yco~X+offset(log(rep(T,2*n))), family=quasipoisson()))
        #Roberts formula for s.e.
        rl$co_quasipoisson.se[i] <- co_quasipoisson.summary$coeffi[2,2]*exp(co_quasipoisson.summary$coeffi[2,1])
        rl$co_quasipoisson.se2[i] <- co_quasipoisson.summary$coeffi[2,2]
        rl$co_quasipoisson.rratio[i] <- exp(co_quasipoisson.summary$coeffi[2,1])   
        rl$co_quasipoisson.pval[i] <- co_quasipoisson.summary$coef[2,4] 
      }
    
      if (('all' %in% meth)|('nbd' %in% meth)) {
        co_nbd.summary <- summary(glm.nb(Yco~X+offset(log(rep(T,2*n)))))
        rl$co_nbd.warn[i] <- (!is.null(co_nbd.summary$th.warn))
        #Roberts formula for s.e.
        rl$co_nbd.se[i] <- co_nbd.summary$coeffi[2,2]*exp(co_nbd.summary$coeffi[2,1])
        rl$co_nbd.se2[i] <- co_nbd.summary$coeffi[2,2]
        rl$co_nbd.rratio[i] <- exp(co_nbd.summary$coeffi[2,1])   
        rl$co_nbd.pval[i] <- co_nbd.summary$coef[2,4] 
        rl$co_nbd.theta[i] <- co_nbd.summary$theta
        rl$mr.co.p[i]=exp(co_nbd.summary$coeffi[1,1])
        rl$mr.co.t[i]=exp(co_nbd.summary$coeffi[1,1] + co_nbd.summary$coeffi[2,1])
      } 
    }
    
    if (('all' %in% meth2)|('dl' %in% meth2)) {
      if (('all' %in% meth)|('poisson' %in% meth)) {    
       # no imputation (direct likelihood approach)
        dl_poisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=poisson()))
    
        #Roberts formula for s.e.
        rl$dl_poisson.se[i] <- dl_poisson.summary$coeffi[2,2]*exp(dl_poisson.summary$coeffi[2,1])
        rl$dl_poisson.se2[i] <- dl_poisson.summary$coeffi[2,2]
        rl$dl_poisson.rratio[i] <- exp(dl_poisson.summary$coeffi[2,1])   #Rate ratio
        rl$dl_poisson.pval[i] <- dl_poisson.summary$coef[2,4] #p-value
      }
    
      if (('all' %in% meth)|('quasipoisson' %in% meth)) {
        dl_quasipoisson.summary <- summary(glm(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0])),family=quasipoisson()))
      
        #Roberts formula for s.e.
        rl$dl_quasipoisson.se[i] <- dl_quasipoisson.summary$coeffi[2,2]*exp(dl_quasipoisson.summary$coeffi[2,1])
        rl$dl_quasipoisson.se2[i] <- dl_quasipoisson.summary$coeffi[2,2]
        rl$dl_quasipoisson.rratio[i] <- exp(dl_quasipoisson.summary$coeffi[2,1])   
        rl$dl_quasipoisson.pval[i] <- dl_quasipoisson.summary$coef[2,4] 
      }
    
      if (('all' %in% meth)|('nbd' %in% meth)) {    
        dl_nbd.summary <- summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))
        rl$dl_nbd.warn[i] <- (!is.null(dl_nbd.summary$th.warn))
          
        #Roberts formula for s.e.
        rl$dl_nbd.se[i] <- dl_nbd.summary$coeffi[2,2]*exp(dl_nbd.summary$coeffi[2,1])
        rl$dl_nbd.se2[i] <- dl_nbd.summary$coeffi[2,2]
        rl$dl_nbd.rratio[i] <- exp(dl_nbd.summary$coeffi[2,1])   
        rl$dl_nbd.pval[i] <- dl_nbd.summary$coef[2,4] 
       
        rl$dl_nbd.theta[i] <- dl_nbd.summary$theta
        rl$mr.dl.p[i]=exp(dl_nbd.summary$coeffi[1,1])
        rl$mr.dl.t[i]=exp(dl_nbd.summary$coeffi[1,1] + dl_nbd.summary$coeffi[2,1])
      }
    }
    
    if (('all' %in% meth2)|('j2r' %in% meth2)) {
      j2r <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 1, N.mi, equal_dispersion=onepertrt,imp_choice2=meth)
      if (('all' %in% meth)|('poisson' %in% meth)) {
        # jump to reference
        rl$j2r_poisson.se[i] <- (j2r$mi_poisson.se2)*exp(j2r$mi_poisson.trt2)
        rl$j2r_poisson.se2[i] <- j2r$mi_poisson.se2
        rl$j2r_poisson.rratio[i] <- exp(j2r$mi_poisson.trt2)
        rl$j2r_poisson.df[i] <- j2r$mi_poisson.df2
        rl$j2r_poisson.adf[i] <- j2r$mi_poisson.adf2
        rl$j2r_poisson.pval[i] <- 2*(1-pt(abs(log(rl$j2r_poisson.rratio[i])/rl$j2r_poisson.se2[i]),df=rl$j2r_poisson.df[i]))
        rl$j2r_poisson.apval[i] <- 2*(1-pt(abs(log(rl$j2r_poisson.rratio[i])/rl$j2r_poisson.se2[i]),df=rl$j2r_poisson.adf[i]))
      }
    
      if (('all' %in% meth)|('quasipoisson' %in% meth)) {
        rl$j2r_quasipoisson.se[i] <- (j2r$mi_quasipoisson.se2)*exp(j2r$mi_quasipoisson.trt2)
        rl$j2r_quasipoisson.se2[i] <- j2r$mi_quasipoisson.se2
        rl$j2r_quasipoisson.rratio[i] <- exp(j2r$mi_quasipoisson.trt2)
        rl$j2r_quasipoisson.df[i] <- j2r$mi_quasipoisson.df2
        rl$j2r_quasipoisson.adf[i] <- j2r$mi_quasipoisson.adf2
        rl$j2r_quasipoisson.pval[i] <- 2*(1-pt(abs(log(rl$j2r_quasipoisson.rratio[i])/rl$j2r_quasipoisson.se2[i]),df=rl$j2r_quasipoisson.df[i]))
        rl$j2r_quasipoisson.apval[i] <- 2*(1-pt(abs(log(rl$j2r_quasipoisson.rratio[i])/rl$j2r_quasipoisson.se2[i]),df=rl$j2r_quasipoisson.adf[i]))
      }
    
      if (('all' %in% meth)|('nbd' %in% meth)) {
        rl$j2r_nbd.se[i] <- (j2r$mi_nbd.se2)*exp(j2r$mi_nbd.trt2)
        rl$j2r_nbd.se2[i] <- j2r$mi_nbd.se2
        rl$j2r_nbd.rratio[i] <- exp(j2r$mi_nbd.trt2)
        rl$j2r_nbd.df[i] <- j2r$mi_nbd.df2
        rl$j2r_nbd.adf[i] <- j2r$mi_nbd.adf2
        rl$j2r_nbd.pval[i] <- 2*(1-pt(abs(log(rl$j2r_nbd.rratio[i])/rl$j2r_nbd.se2[i]),df=rl$j2r_nbd.df[i]))
        rl$j2r_nbd.apval[i] <- 2*(1-pt(abs(log(rl$j2r_nbd.rratio[i])/rl$j2r_nbd.se2[i]),df=rl$j2r_nbd.adf[i]))
        rl$j2r_nbd.th[i] <- j2r$mi_nbd.th 
        rl$mr.j2r.p[i]=j2r$mr.p
        rl$events.j2r.p[i]=j2r$ev.p
    
        rl$mr.j2r.t[i]=j2r$mr.t
        rl$events.j2r.t[i]=j2r$ev.t
        rl$p.rate[i]=j2r$p.rate
        rl$t.rate[i]=j2r$t.rate
      }
    }  
    # efficacy
    if (('all' %in% meth2)|('efy' %in% meth2)) {
      efy <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 2, N.mi, equal_dispersion=onepertrt,imp_choice2=meth)
      if (('all' %in% meth)|('poisson' %in% meth)) {
        rl$efy_poisson.se[i] <- (efy$mi_poisson.se2)*exp(efy$mi_poisson.trt2)
        rl$efy_poisson.se2[i] <- efy$mi_poisson.se2
        rl$efy_poisson.rratio[i] <- exp(efy$mi_poisson.trt2)
        rl$efy_poisson.df[i] <- efy$mi_poisson.df2
        rl$efy_poisson.adf[i] <- efy$mi_poisson.adf2
        rl$efy_poisson.pval[i] <- 2*(1-pt(abs(log(rl$efy_poisson.rratio[i])/rl$efy_poisson.se2[i]),df=rl$efy_poisson.df[i]))
        rl$efy_poisson.apval[i] <- 2*(1-pt(abs(log(rl$efy_poisson.rratio[i])/rl$efy_poisson.se2[i]),df=rl$efy_poisson.adf[i]))
      }
      if (('all' %in% meth)|('quasipoisson' %in% meth)) {
        rl$efy_quasipoisson.se[i] <- (efy$mi_quasipoisson.se2)*exp(efy$mi_quasipoisson.trt2)
        rl$efy_quasipoisson.se2[i] <- efy$mi_quasipoisson.se2
        rl$efy_quasipoisson.rratio[i] <- exp(efy$mi_quasipoisson.trt2)
        rl$efy_quasipoisson.df[i] <- efy$mi_quasipoisson.df2
        rl$efy_quasipoisson.adf[i] <- efy$mi_quasipoisson.adf2
        rl$efy_quasipoisson.pval[i] <- 2*(1-pt(abs(log(rl$efy_quasipoisson.rratio[i])/rl$efy_quasipoisson.se2[i]),df=rl$efy_quasipoisson.df[i]))
        rl$efy_quasipoisson.apval[i] <- 2*(1-pt(abs(log(rl$efy_quasipoisson.rratio[i])/rl$efy_quasipoisson.se2[i]),df=rl$efy_quasipoisson.adf[i]))
      }
      if (('all' %in% meth)|('nbd' %in% meth)) {
        rl$efy_nbd.se[i] <- (efy$mi_nbd.se2)*exp(efy$mi_nbd.trt2)
        rl$efy_nbd.se2[i] <- efy$mi_nbd.se2
        rl$efy_nbd.rratio[i] <- exp(efy$mi_nbd.trt2)
        rl$efy_nbd.df[i] <- efy$mi_nbd.df2
        rl$efy_nbd.adf[i] <- efy$mi_nbd.adf2
        rl$efy_nbd.pval[i] <- 2*(1-pt(abs(log(rl$efy_nbd.rratio[i])/rl$efy_nbd.se2[i]),df=rl$efy_nbd.df[i]))
        rl$efy_nbd.apval[i] <- 2*(1-pt(abs(log(rl$efy_nbd.rratio[i])/rl$efy_nbd.se2[i]),df=rl$efy_nbd.adf[i]))
      
        rl$mr.ef.p[i]=efy$mr.p
        rl$events.ef.p[i]=efy$ev.p
        rl$mr.ef.t[i]=efy$mr.t
        rl$events.ef.t[i]=efy$ev.t
        #if (!is.null(efy$outliers)) outliers.ef=rbind(outliers.ef, cbind(efy$outliers, i))
    
        rl$p.rate[i]=efy$p.rate
        rl$t.rate[i]=efy$t.rate
      }
    }
    
    
    if (('all' %in% meth2)|('tr' %in% meth2)) {
      tr <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 3, N.mi, equal_dispersion=onepertrt,imp_choice2=meth, delta=delta)
      if (('all' %in% meth)|('poisson' %in% meth)) {
        # jump to reference
        rl$tr_poisson.se[i] <- (tr$mi_poisson.se2)*exp(tr$mi_poisson.trt2)
        rl$tr_poisson.se2[i] <- tr$mi_poisson.se2
        rl$tr_poisson.rratio[i] <- exp(tr$mi_poisson.trt2)
        rl$tr_poisson.df[i] <- tr$mi_poisson.df2
        rl$tr_poisson.adf[i] <- tr$mi_poisson.adf2
        rl$tr_poisson.pval[i] <- 2*(1-pt(abs(log(rl$tr_poisson.rratio[i])/rl$tr_poisson.se2[i]),df=rl$tr_poisson.df[i]))
        rl$tr_poisson.apval[i] <- 2*(1-pt(abs(log(rl$tr_poisson.rratio[i])/rl$tr_poisson.se2[i]),df=rl$tr_poisson.adf[i]))
      }
    
      if (('all' %in% meth)|('quasipoisson' %in% meth)) {
        rl$tr_quasipoisson.se[i] <- (tr$mi_quasipoisson.se2)*exp(tr$mi_quasipoisson.trt2)
        rl$tr_quasipoisson.se2[i] <- tr$mi_quasipoisson.se2
        rl$tr_quasipoisson.rratio[i] <- exp(tr$mi_quasipoisson.trt2)
        rl$tr_quasipoisson.df[i] <- tr$mi_quasipoisson.df2
        rl$tr_quasipoisson.adf[i] <- tr$mi_quasipoisson.adf2
        rl$tr_quasipoisson.pval[i] <- 2*(1-pt(abs(log(rl$tr_quasipoisson.rratio[i])/rl$tr_quasipoisson.se2[i]),df=rl$tr_quasipoisson.df[i]))
        rl$tr_quasipoisson.apval[i] <- 2*(1-pt(abs(log(rl$tr_quasipoisson.rratio[i])/rl$tr_quasipoisson.se2[i]),df=rl$tr_quasipoisson.adf[i]))
      }
    
      if (('all' %in% meth)|('nbd' %in% meth)) {
        rl$tr_nbd.se[i] <- (tr$mi_nbd.se2)*exp(tr$mi_nbd.trt2)
        rl$tr_nbd.se2[i] <- tr$mi_nbd.se2
        rl$tr_nbd.rratio[i] <- exp(tr$mi_nbd.trt2)
        rl$tr_nbd.df[i] <- tr$mi_nbd.df2
        rl$tr_nbd.adf[i] <- tr$mi_nbd.adf2
        rl$tr_nbd.pval[i] <- 2*(1-pt(abs(log(rl$tr_nbd.rratio[i])/rl$tr_nbd.se2[i]),df=rl$tr_nbd.df[i]))
        rl$tr_nbd.apval[i] <- 2*(1-pt(abs(log(rl$tr_nbd.rratio[i])/rl$tr_nbd.se2[i]),df=rl$tr_nbd.adf[i]))
        
        rl$mr.tr.p[i]=tr$mr.p
        rl$events.tr.p[i]=tr$ev.p
    
        rl$mr.tr.t[i]=tr$mr.t
        rl$events.tr.t[i]=tr$ev.t
        
        rl$p.rate[i]=tr$p.rate
        rl$t.rate[i]=tr$t.rate
      }
    }
    
    if (('all' %in% meth2)|('delta' %in% meth2)) {
      del <- cim(T, n, nb_complete.placebo, nb_incomplete.placebo, nb_complete.treatment, nb_incomplete.treatment, 4, N.mi, equal_dispersion=onepertrt,imp_choice2=meth, delta=delta)
      if (('all' %in% meth)|('poisson' %in% meth)) {
        # jump to reference
        rl$delta_poisson.se[i] <- (del$mi_poisson.se2)*exp(del$mi_poisson.trt2)
        rl$delta_poisson.se2[i] <- del$mi_poisson.se2
        rl$delta_poisson.rratio[i] <- exp(del$mi_poisson.trt2)
        rl$delta_poisson.df[i] <- del$mi_poisson.df2
        rl$delta_poisson.adf[i] <- del$mi_poisson.adf2
        rl$delta_poisson.pval[i] <- 2*(1-pt(abs(log(rl$delta_poisson.rratio[i])/rl$delta_poisson.se2[i]),df=rl$delta_poisson.df[i]))
        rl$delta_poisson.apval[i] <- 2*(1-pt(abs(log(rl$delta_poisson.rratio[i])/rl$delta_poisson.se2[i]),df=rl$delta_poisson.adf[i]))
      }
    
      if (('all' %in% meth)|('quasipoisson' %in% meth)) {
        rl$delta_quasipoisson.se[i] <- (del$mi_quasipoisson.se2)*exp(del$mi_quasipoisson.trt2)
        rl$delta_quasipoisson.se2[i] <- del$mi_quasipoisson.se2
        rl$delta_quasipoisson.rratio[i] <- exp(del$mi_quasipoisson.trt2)
        rl$delta_quasipoisson.df[i] <- del$mi_quasipoisson.df2
        rl$delta_quasipoisson.adf[i] <- del$mi_quasipoisson.adf2
        rl$delta_quasipoisson.pval[i] <- 2*(1-pt(abs(log(rl$delta_quasipoisson.rratio[i])/rl$delta_quasipoisson.se2[i]),df=rl$delta_quasipoisson.df[i]))
        rl$delta_quasipoisson.apval[i] <- 2*(1-pt(abs(log(rl$delta_quasipoisson.rratio[i])/rl$delta_quasipoisson.se2[i]),df=rl$delta_quasipoisson.adf[i]))
      }
    
      if (('all' %in% meth)|('nbd' %in% meth)) {
        rl$delta_nbd.se[i] <- (del$mi_nbd.se2)*exp(del$mi_nbd.trt2)
        rl$delta_nbd.se2[i] <- del$mi_nbd.se2
        rl$delta_nbd.rratio[i] <- exp(del$mi_nbd.trt2)
        rl$delta_nbd.df[i] <- del$mi_nbd.df2
        rl$delta_nbd.adf[i] <- del$mi_nbd.adf2
        rl$delta_nbd.pval[i] <- 2*(1-pt(abs(log(rl$delta_nbd.rratio[i])/rl$delta_nbd.se2[i]),df=rl$delta_nbd.df[i]))
        rl$delta_nbd.apval[i] <- 2*(1-pt(abs(log(rl$delta_nbd.rratio[i])/rl$delta_nbd.se2[i]),df=rl$delta_nbd.adf[i]))
    
        rl$delta_nbd.th[i] <- del$mi_nbd.th
        rl$mr.delta.p[i]=del$mr.p
        rl$events.delta.p[i]=del$ev.p
    
        rl$mr.delta.t[i]=del$mr.t
        rl$events.delta.t[i]=del$ev.t
        rl$p.rate[i]=del$p.rate
        rl$t.rate[i]=del$t.rate
      }
    }
    
    rl$err[i] <- 0
    }, warning = function(w){return(1)})
    return(list(rl=rl))
}



comb <- function(x, ...) {  
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}


sim_mar2_mc<-function(n, rate.placebo, rr ,T ,dispersion.placebo, dispersion.treatment,rates,
                      drop_rate_variance,N.rep,N.mi, alpha, onepertrt, datname, meth, meth2, 
                      delta=2){
  startdf=setdf(1)  
     
  outliers.ef=c()
  outliers.j2r=c()
  
  rate.treatment <- (1-rr)*rate.placebo
  X <- c(rep(0,n), rep(1,n))
  cdat <- genNmar(N.rep, n=n,rate.placebo=rate.placebo,rr=rr,T=T,
                      dispersion.placebo=dispersion.placebo, dispersion.treatment=dispersion.treatment,
                      rates=rates,drop_rate_variance=drop_rate_variance)
    
  res <- foreach(i=1:N.rep, .combine=comb, 
                     .export=c('gc','nb_complete','dropout_time4','nb_incomplete',
                               'subset','glm.nb','cim','makeDataFrame') ) %dopar% 
             {
              gc(cdat,j=i, T,n, N.mi, onepertrt, X, startdf, meth, meth2, delta=delta) 
              #rate.placebo,dispersion.placebo, rate.treatment,dispersion.treatment, rates, drop_rate_variance,
             }
  attach(res$rl)
  # excludes cases with warnings
  co_powerpois<-sum(subset(co_poisson.pval,err==0)<alpha)/length(subset(co_poisson.pval,err==0))
  co_powerqpois<-sum(subset(co_quasipoisson.pval,err==0)<alpha)/length(subset(co_quasipoisson.pval,err==0))
  co_powernb<-sum(subset(co_nbd.pval,err==0)<alpha)/length(subset(co_nbd.pval,err==0))

  dl_powerpois <-sum(subset(dl_poisson.pval,     err==0)<alpha)/length(subset(dl_poisson.pval,err==0)) 
  dl_powerqpois<-sum(subset(dl_quasipoisson.pval,err==0)<alpha)/length(subset(dl_quasipoisson.pval,err==0)) 
  dl_powernb   <-sum(subset(dl_nbd.pval,         err==0)<alpha)/length(subset(dl_nbd.pval,err==0)) 
  
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
  
  tr_powerpois<-sum(subset(tr_poisson.pval,err==0)<alpha)/length(subset(tr_poisson.pval,err==0))
  tr_powerqpois<-sum(subset(tr_quasipoisson.pval,err==0)<alpha)/length(subset(tr_quasipoisson.pval,err==0))
  tr_powernb<-sum(subset(tr_nbd.pval,err==0)<alpha)/length(subset(tr_nbd.pval,err==0))
  
  tr_powerpoisa<-sum(subset(tr_poisson.apval,err==0)<alpha)/length(subset(tr_poisson.apval,err==0))
  tr_powerqpoisa<-sum(subset(tr_quasipoisson.apval,err==0)<alpha)/length(subset(tr_quasipoisson.apval,err==0))
  tr_powernba<-sum(subset(tr_nbd.apval,err==0)<alpha)/length(subset(tr_nbd.apval,err==0))
  
  del_powerpois<-sum(subset(delta_poisson.pval,err==0)<alpha)/length(subset(delta_poisson.pval,err==0))
  del_powerqpois<-sum(subset(delta_quasipoisson.pval,err==0)<alpha)/length(subset(delta_quasipoisson.pval,err==0))
  del_powernb<-sum(subset(delta_nbd.pval,err==0)<alpha)/length(subset(delta_nbd.pval,err==0))
  
  del_powerpoisa<-sum(subset(delta_poisson.apval,err==0)<alpha)/length(subset(delta_poisson.apval,err==0))
  del_powerqpoisa<-sum(subset(delta_quasipoisson.apval,err==0)<alpha)/length(subset(delta_quasipoisson.apval,err==0))
  del_powernba<-sum(subset(delta_nbd.apval,err==0)<alpha)/length(subset(delta_nbd.apval,err==0))
 
  
  power <- as.data.frame(cbind(co_powerpois,co_powerqpois,co_powernb,
                               dl_powerpois,dl_powerqpois,dl_powernb,
                               j2r_powerpois,j2r_powerqpois,j2r_powernb,j2r_powerpoisa,j2r_powerqpoisa,j2r_powernba,
                               efy_powerpois,efy_powerqpois,efy_powernb,efy_powerpoisa,efy_powerqpoisa,efy_powernba,
                               tr_powerpois,tr_powerqpois,tr_powernb,tr_powerpoisa,tr_powerqpoisa,tr_powernba,
                               del_powerpois,del_powerqpois,del_powernb,del_powerpoisa,del_powerqpoisa,del_powernba))
 
  rratio <- as.data.frame(cbind(co_poisson.rratio,co_quasipoisson.rratio,co_nbd.rratio,
                                dl_poisson.rratio,dl_quasipoisson.rratio,dl_nbd.rratio,
                                j2r_poisson.rratio,j2r_quasipoisson.rratio,j2r_nbd.rratio,
                                efy_poisson.rratio,efy_quasipoisson.rratio,efy_nbd.rratio,
                                tr_poisson.rratio,tr_quasipoisson.rratio,tr_nbd.rratio,
                                delta_poisson.rratio,delta_quasipoisson.rratio,delta_nbd.rratio))
  
  se <- as.data.frame(cbind(co_poisson.se,co_quasipoisson.se,co_nbd.se,
                            dl_poisson.se,dl_quasipoisson.se,dl_nbd.se,
                            j2r_poisson.se,j2r_quasipoisson.se,j2r_nbd.se,
                            efy_poisson.se,efy_quasipoisson.se,efy_nbd.se,
                            tr_poisson.se,tr_quasipoisson.se,tr_nbd.se,
                            delta_poisson.se,delta_quasipoisson.se,delta_nbd.se,
                            co_poisson.se2,co_quasipoisson.se2,co_nbd.se2,
                            dl_poisson.se2,dl_quasipoisson.se2,dl_nbd.se2,
                            j2r_poisson.se2,j2r_quasipoisson.se2,j2r_nbd.se2,
                            efy_poisson.se2,efy_quasipoisson.se2,efy_nbd.se2,
                            tr_poisson.se2,tr_quasipoisson.se2,tr_nbd.se2,
                            delta_poisson.se2,delta_quasipoisson.se2,delta_nbd.se2))
  
  pval <- as.data.frame(cbind(co_poisson.pval,co_quasipoisson.pval,co_nbd.pval,
                            dl_poisson.pval,dl_quasipoisson.pval,dl_nbd.pval,
                            j2r_poisson.pval,j2r_quasipoisson.pval,j2r_nbd.pval,
                            efy_poisson.pval,efy_quasipoisson.pval,efy_nbd.pval,
                            tr_poisson.pval,tr_quasipoisson.pval,tr_nbd.pval,
                            delta_poisson.pval,delta_quasipoisson.pval,delta_nbd.pval,
                            co_poisson.apval,co_quasipoisson.apval,co_nbd.apval,
                            dl_poisson.apval,dl_quasipoisson.apval,dl_nbd.apval,
                            j2r_poisson.apval,j2r_quasipoisson.apval,j2r_nbd.apval,
                            efy_poisson.apval,efy_quasipoisson.apval,efy_nbd.apval,
                            tr_poisson.apval,tr_quasipoisson.apval,tr_nbd.apval,
                            delta_poisson.apval,delta_quasipoisson.apval,delta_nbd.apval))
  
  mi_df <- as.data.frame(cbind(j2r_poisson.df,j2r_poisson.adf,j2r_quasipoisson.df,j2r_quasipoisson.adf,j2r_nbd.df,j2r_nbd.adf,
                               efy_poisson.df,efy_poisson.adf,efy_quasipoisson.df,efy_quasipoisson.adf,efy_nbd.df,efy_nbd.adf,
                               tr_poisson.df,tr_poisson.adf,tr_quasipoisson.df,tr_quasipoisson.adf,tr_nbd.df,tr_nbd.adf,
                               delta_poisson.df,delta_poisson.adf,delta_quasipoisson.df,delta_quasipoisson.adf,delta_nbd.df,delta_nbd.adf))
  
  warnalltype <- sum(err==1)
  
  summdrop.p <- summary(dropout.p)
  summdrop.t <- summary(dropout.t)
  dropout <- as.data.frame(cbind(summdrop.p,summdrop.t))
  
  censortime <- as.data.frame(cbind(P=summary(censor.time.p),T=summary(censor.time.t)))
  censortimeT <- as.data.frame(cbind(P=summary(censor.timeT.p),T=summary(censor.timeT.t)))
  
  ERgt75p <- cbind(P=ERgt75p.p,T=ERgt75p.t,RR=ERgt75p.t/ERgt75p.p)
  ERle75p <- cbind(P=ERle75p.p,T=ERle75p.t,RR=ERle75p.t/ERle75p.p)
  
  write.csv(cdat$dat,datname)
  
  result <- list(err=err,power=power, rratio=rratio,se=se,mi_df=mi_df,pval=pval,dropout=dropout,censortime=censortime,censortimeT=censortimeT,ERgt75p=ERgt75p,ERle75p=ERle75p, 
                 mr.dl.p=mr.dl.p, mr.ef.p=mr.ef.p, mr.co.p=mr.co.p, mr.j2r.p=mr.j2r.p,mr.tr.p=mr.tr.p, mr.delta.p=mr.delta.p, 
                 events.dl.p=events.dl.p, events.ef.p=events.ef.p, events.j2r.p=events.j2r.p, events.co.p=events.co.p, events.tr.p=events.tr.p, events.delta.p=events.delta.p,
                 mr.dl.t=mr.dl.t, mr.ef.t=mr.ef.t, mr.co.t=mr.co.t, mr.j2r.t=mr.j2r.t, mr.tr.t=mr.tr.t, mr.delta.t=mr.delta.t, 
                 events.dl.t=events.dl.t, events.ef.t=events.ef.t, events.j2r.t=events.j2r.t, events.co.t=events.co.t, events.tr.t=events.tr.t, events.delta.t=events.delta.t, 
                 warnalltype=warnalltype, j2r_nbd.th  = j2r_nbd.th, delta_nbd.th  = delta_nbd.th)
  detach(res$rl)
  return(result)
}

