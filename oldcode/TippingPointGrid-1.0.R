#### Graphical sensitivity analysis for recurrent events endpoint.
#### Similar to what's described in
#### Sally Hollis (2002) - A graphical sensitivity analysis for clinical trials with non-ignorable missing binary outcome 
#### Victoria Liublinska and Donald B. Rubin (2013) - Sensitivity analysis for a partially missing binary outcome in a 
####                                                  two-arm randomized clinical trial
####
#### Mattis Gottlow, 2015-06-18
require(ggplot2)
require(MASS)
path="\\\\emea.astrazeneca.net\\sweden\\Molndal\\Users 05\\KDSS481\\Documents\\Training\\Missing data\\INR"
source(paste(path,"/Code/NBgenfun.R", sep=""))
source(paste(path,"/Code/power_functionsv2.R", sep=""))
source(paste(path,"/Code/cim_v1.5.R", sep=""))

####### Controlled Imputation function ########
 cimgrid=function(d1,d2,n.mi) {
  tr=cim(1, n, nb_complete.placebo, nb_incomplete.placebo, 
              nb_complete.treatment, nb_incomplete.treatment, 
              imp_choice=4, n.mi, equal_dispersion=1, imp_choice2="nbd", delta=c(d1,d2))  
  se2 <- tr$mi_nbd.se2
  rratio <- exp(tr$mi_nbd.trt2)
  adf <- tr$mi_nbd.adf2
  apval<- 2*(1-pt(abs(log(rratio)/se2),df=adf))
  return(data.frame(d1,d2,rratio, tr$ev.p, tr$ev.t, tr$mr.p, tr$mr.t,tr$p.rate, tr$t.rate,se2, pval=apval))
 }
    
##############  Normal approximation function #################
nbtst = function(alpha,c0,c1,k,n){
  est=log(c1/c0)
  V=(1/c0 + 1/c1 + 2*k/n)
  pval=2*(1-pnorm(abs(est/sqrt(V))))
  return(data.frame(pval=pval,effect=exp(est),V=V))
}

#Generate data set and create variables needed
n=250
cdati=genNmar(1,n,0.9,0.4,1,1,1,seq(0.05,20,1),0.05)
cdi.p=cdati$dat[cdati$dat$TRTn==0,]
cdi.t=cdati$dat[cdati$dat$TRTn==1,]
nb_inc.p <- cdi.p$COUNT
nb_inc.t <- cdi.t$COUNT
do.p <- cdi.p$STUDY_TIME
do.t <- cdi.t$STUDY_TIME
nb_complete.placebo=cdi.p$cc
nb_complete.treatment=cdi.t$cc  
Y <- cdati$dat$COUNT
X <- c(rep(0,n),rep(1,n))
censor.time <- cdati$dat$STUDY_TIME

nbsum=summary(glm.nb(Y[censor.time>0]~X[censor.time>0]+offset(log(censor.time[censor.time>0]))))  
theta=nbsum$theta

nb_incomplete.placebo=list(nb_inc.p, do.p, cdati$times.p[(j-1)*n + 1:n])
nb_incomplete.treatment=list(nb_inc.t, do.t, cdati$times.t[(j-1)*n + 1:n])
    
py.p=sum(do.p)
py.t=sum(do.t)
ev.p=sum(nb_inc.p)
ev.t=sum(nb_inc.t)
r.p=ev.p/py.p
r.t=ev.t/py.t

####### Convert deltas to number of events ####### 
cg=expand.grid(d1=(8:20)/10,d2=(10:22)/10)
cg$imp0 = round(cg$d1*(r.p)*(n-py.p)/theta)
cg$imp1 = round(cg$d2*(r.t)*(n-py.t)/theta)

####### Normal approximation data #######
gridres1=data.frame(d1=cg$d1,d2=cg$d2,cg$imp0,cg$imp1, rr=ev.t/ev.p,
                   nbtst(0.05,ev.p+cg$imp0,ev.t+cg$imp1, 1/theta, n))
### Plot ###
p <- ggplot(gridres1, aes(d1,d2, fill=pval*(pval>0.05))) + 
  geom_tile(colour = "white") + #geom_tile(aes(alpha = pval)
  scale_fill_gradient(low = "green", high = "red", guide=FALSE)
p = p + ylab("Delta of Active arm") + 
        xlab("Delta of Placebo arm")
p 
 
 
####### COntrolled Imputation data #######
gridres2=data.frame()
for (i in 1:length(cg$d1)) {
     print(c(i,length(cg$d1)))
     gridresi=cimgrid(cg$d1[i],cg$d2[i],30)
     gridres2=rbind(gridres2,gridresi)
}
### Plot ###
p <- ggplot(gridres2, aes(d1,d2, fill=pval*(pval>0.05))) + 
  geom_tile(colour = "white") + #geom_tile(aes(alpha = pval)
  scale_fill_gradient(low = "green", high = "red", guide=FALSE)
p = p + ylab("Delta of Active arm") + 
        xlab("Delta of Placebo arm")
p
  
### Plot with bothe methods ###
gr1=gridres1[,c(1,2,6)]  
gr2=gridres2[,c(1,2,11)]
gr1$method="Normal Approximation"
gr2$method="Controlled Imputation"
gr=rbind(gr1,gr2)  
    
p2 <- ggplot(gr, aes(d1,d2, fill=pval*(pval>0.05))) + 
  geom_tile(colour = "white") + #geom_tile(aes(alpha = pval)
  scale_fill_gradient(low = "green", high = "red", guide=FALSE)
p2 = p2 + facet_wrap(~method)    
p2 = p2 + ylab("Delta of Active arm") + 
        xlab("Delta of Placebo arm")
p2
    
    