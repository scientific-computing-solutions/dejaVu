
path="\\\\emea.astrazeneca.net\\sweden\\Molndal\\Users 05\\KDSS481\\Documents\\Training\\Missing data\\INR"
source(paste(path,"/Code/NBgenfun.R", sep=""))
source(paste(path,"/Code/power_functionsv2.R", sep=""))
source(paste(path,"/Code/cim_v1.5.R", sep=""))

library(doSNOW)
library(foreach)
cl <- makeCluster(3)
registerDoSNOW(cl)
nr=2000
library(ggplot2)
getn = function(dat,k) return(dat[!is.na(dat[,k]),])

# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
drates=seq(0.001,5,0.1)
sasname=paste(path,"/Data/PSIposter1.csv",sep="")
lstname=paste(path,"/Data/PSIposter1.rds",sep="")
ptime <- system.time({
poster1 <- sim_mar2_mc(n=228,rate.placebo=0.88,rr=0.4,T=1,
                      dispersion.placebo=0.9, dispersion.treatment=0.9,
                      rates=drates,drop_rate_variance=0.05,
                      N.rep=nr, N.mi=10, alpha=0.04, 1, sasname, meth=c("nbd"),
			    meth2=c("complete","dl","j2r","tr","efy","delta"),delta=c(1,1.4))
})[3]
ptime

poster1$power
poster1=within(poster1, se$co_nbd.se_adj<-se$co_nbd.se/sqrt((py.p+py.t)/2))
saveRDS(poster1, file=lstname)


# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
drates=seq(0.001,10,0.25)
sasname=paste(path,"/Data/PSIposter2.csv",sep="")
lstname=paste(path,"/Data/PSIposter2.rds",sep="")
ptime <- system.time({
poster2 <- sim_mar2_mc(n=228,rate.placebo=0.88,rr=0.4,T=1,
                             dispersion.placebo=0.9,dispersion.treatment=0.9,
                             rates=drates,drop_rate_variance=0.05,
                             N.rep=nr,N.mi=10,alpha=0.04,1,sasname,meth=c("nbd"), 
				     meth2=c("complete","dl","j2r","tr","efy","delta"),delta=c(1,1.4))
})[3]
ptime
poster2=within(poster2, se$co_nbd.se_adj<-se$co_nbd.se/sqrt((py.p+py.t)/2))
saveRDS(poster2, file=lstname)


# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
drates=seq(0.001,20,0.5)
sasname=paste(path,"/Data/PSIposter3.csv",sep="")
lstname=paste(path,"/Data/PSIposter3.rds",sep="")
ptime <- system.time({
poster3 <- sim_mar2_mc(n=228,rate.placebo=0.88,rr=0.4,T=1,
                             dispersion.placebo=0.9,dispersion.treatment=0.9,
                             rates=drates,drop_rate_variance=0.05,
                             N.rep=nr,N.mi=10,alpha=0.04,1,sasname,meth=c("nbd"), 
				     meth2=c("complete","dl","j2r","tr","efy","delta"),delta=c(1,1.4))
})[3]
ptime
poster3=within(poster3, se$co_nbd.se_adj<-se$co_nbd.se/sqrt((py.p+py.t)/2))
saveRDS(poster3, file=lstname)

# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
drates=seq(0.001,30,1)
sasname=paste(path,"/Data/PSIposter4.csv",sep="")
lstname=paste(path,"/Data/PSIposter4.rds",sep="")
ptime <- system.time({
poster4 <- sim_mar2_mc(n=228,rate.placebo=0.88,rr=0.4,T=1,
                             dispersion.placebo=0.9,dispersion.treatment=0.9,
                             rates=drates,drop_rate_variance=0.05,
                             N.rep=nr,N.mi=10,alpha=0.04,1,sasname,meth=c("nbd"), 
				     meth2=c("complete","dl","j2r","tr","efy","delta"),delta=c(1,1.4))
})[3]
ptime
poster4=within(poster4, se$co_nbd.se_adj<-se$co_nbd.se/sqrt((py.p+py.t)/2))
saveRDS(poster4, file=lstname)


# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
drates=seq(0.001,40,1.5)
sasname=paste(path,"/Data/PSIposter5.csv",sep="")
lstname=paste(path,"/Data/PSIposter5.rds",sep="")
ptime <- system.time({
poster5 <- sim_mar2_mc(n=228,rate.placebo=0.88,rr=0.4,T=1,
                             dispersion.placebo=0.9,dispersion.treatment=0.9,
                             rates=drates,drop_rate_variance=0.05,
                             N.rep=nr,N.mi=10,alpha=0.04,1,sasname,meth=c("nbd"), 
				     meth2=c("complete","dl","j2r","tr","efy","delta"),delta=c(1,1.4))
})[3]
ptime
poster5=within(poster5, se$co_nbd.se_adj<-se$co_nbd.se/sqrt((py.p+py.t)/2))
saveRDS(poster5, file=lstname)


# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
drates=seq(0.001,40,2)
sasname=paste(path,"/Data/PSIposter6.csv",sep="")
lstname=paste(path,"/Data/PSIposter6.rds",sep="")
ptime <- system.time({
poster6 <- sim_mar2_mc(n=228,rate.placebo=0.88,rr=0.4,T=1,
                             dispersion.placebo=0.9,dispersion.treatment=0.9,
                             rates=drates,drop_rate_variance=0.05,
                             N.rep=nr,N.mi=10,alpha=0.04,1,sasname,meth=c("nbd"), 
				     meth2=c("complete","dl","j2r","tr","efy","delta"),delta=c(1,1.4))
})[3]
ptime
poster6=within(poster6, se$co_nbd.se_adj<-se$co_nbd.se/sqrt((py.p+py.t)/2))
saveRDS(poster6, file=lstname)




#tst=readRDS(lstname)
#tst=readRDS("//samba-hpc/raac_scratch/INR/missingData/Data/PSIposter5.rds")

#some checks
mean(poster5$mr.delta.t/poster5$mr.delta.p)
colMeans(poster5$rratio) #stämmer
poster5$censortimeT
#förväntad active rate (lägre än observerad
0.88*0.6*0.9+0.88*0.6*3*0.1
mean(poster5$mr.delta.t)
mean(poster5$mr.delta.p)
mean(poster5$mr.delta.t/poster5$mr.delta.p)
colMeans(poster5$se) #Mycket högre se för delta metod
thv=rbind(c(mean(1/poster1$delta_nbd.th), mean(1/poster1$j2r_nbd.th)),
c(mean(1/poster2$delta_nbd.th), mean(1/poster2$j2r_nbd.th)),
c(mean(1/poster3$delta_nbd.th), mean(1/poster3$j2r_nbd.th)),
c(mean(1/poster4$delta_nbd.th), mean(1/poster4$j2r_nbd.th)),
c(mean(1/poster5$delta_nbd.th), mean(1/poster5$j2r_nbd.th)))#Larger for delata than for j2r
thv[,1]/thv[,2]

#Increase mean -> increase variance
#		   -> increase dispersion 


#approximate se given different thetas?
16.2*sqrt(0.7+1.34*(0.7)^2)/sqrt(0.7+1.03*(0.7)^2)

scens=c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Scenario 5")
meths=c("Complete","DL","J2R","Delta")


do=rbind(100*round(poster1$dropout[4,1:2],digits=3),
	   100*round(poster2$dropout[4,1:2],digits=3),
	   100*round(poster3$dropout[4,1:2],digits=3),
	   100*round(poster4$dropout[4,1:2],digits=3),
	   100*round(poster5$dropout[4,1:2],digits=3)
	   )
doq1=rbind(100*round(poster1$dropout[2,1:2],digits=3),
	   100*round(poster2$dropout[2,1:2],digits=3),
	   100*round(poster3$dropout[2,1:2],digits=3),
	   100*round(poster4$dropout[2,1:2],digits=3),
	   100*round(poster5$dropout[2,1:2],digits=3)
	   )

doq2=rbind(100*round(poster1$dropout[5,1:2],digits=3),
	   100*round(poster2$dropout[5,1:2],digits=3),
	   100*round(poster3$dropout[5,1:2],digits=3),
	   100*round(poster4$dropout[5,1:2],digits=3),
	   100*round(poster5$dropout[5,1:2],digits=3)
	   )

sc=cbind(paste0(do[,1], "% (",doq1[,1], "%-", doq2[,1], "%)"),
	   paste0(do[,2], "% (",doq1[,2], "%-", doq2[,2], "%)"))

rownames(sc)=scens
colnames(sc)=c("Placebo\nmean(Q1-Q3)","Active\nmean(Q1-Q3)")
grid.arrange(tableGrob(sc,core.just="left",gpar.rowtext = gpar(cex = 1, fontface = "bold")))


getres=function(dat,coln,rown) {
	allnm=names(dat)
	nbdnm=allnm[!is.na(dat[1,])]#[!is.na(allnm[!is.na(dat)])]
	ind=which(allnm %in% nbdnm)
	nb=round(dat[,ind],digits=4)
	colnames(nb)=coln
	rownames(nb)=rown
	return(nb)
}

allpow=rbind(poster1$power,poster2$power,poster3$power,poster4$power,poster5$power)
allnm=names(allpow)
nbdnm=paste0(c("co_", "dl_", "j2r_", "del_" ),"powernb")
ind=which(allnm %in% nbdnm)
nb=round(allpow[,ind],digits=4)
colnames(nb)=meths
rownames(nb)=scens
nb

domean=(do[,1]+do[,2])/2
nbg=rbind(data.frame(Power=as.numeric(nb$Complete),Method="Complete",scens=domean),
          data.frame(Power=as.numeric(nb$DL),Method="DL",scens=domean),
	    data.frame(Power=as.numeric(nb$J2R),Method="J2R",scens=domean),
          data.frame(Power=as.numeric(nb$Delta),Method="Delta=2",scens=domean))
colnames(nbg)=c("Power","Method","scens")
nbg
m <- ggplot(nbg, aes(x=scens,y=Power, colour=Method))
m +  geom_line(size=1.5)+ xlab("Dropout Rate")

pow.csp(alpha = 0.025, beta=0.1, mfu = 0.975, lamC = 0.9,  
RRR = c(0.17,0.25,0.3,0.35,0.4), k = 0.95, N_per_arm = 230 )

grid.arrange(main="Type 1 error",
textGrob("Poisson", y=0.2),
tableGrob(pois),
textGrob("QuasiPoisson", y=0.2),
tableGrob(qpois),
textGrob("Negative Binomial", y=0.2),
tableGrob(nb), ncol=1
)


#Error
mse=rbind(colMeans(getn(poster1$se,6)),
	    colMeans(getn(poster2$se,6)),
	    colMeans(getn(poster3$se,6)),
	    colMeans(getn(poster4$se,6)),
	    colMeans(getn(poster5$se,6)))
allnm=colnames(mse)
nbdnm=paste0(c("co_", "dl_", "j2r_", "delta_" ),"nbd.se")
ind=which(allnm %in% nbdnm)
nb.se=round(mse[,ind],digits=4)
colnames(nb.se)=meths
rownames(nb.se)=scens
nb.se

grid.arrange(main="Standard Error (log)",
textGrob("Poisson", y=0.1),
tableGrob(pois.se),
textGrob("QuasiPoisson", y=0.1),
tableGrob(qpois.se),
textGrob("Negative Binomial", y=0.1),
tableGrob(nb.se), ncol=1
)


#AERR
rr=as.data.frame(rbind(colMeans(getn(poster1$rratio,6)),
	   colMeans(getn(poster2$rratio,6)),
	   colMeans(getn(poster3$rratio,6)),
	   colMeans(getn(poster4$rratio,6)),
	   colMeans(getn(poster5$rratio,6))))
allnm=names(rr)
nbdnm=paste0(c("co_", "dl_", "j2r_", "delta_" ),"nbd.rratio")
ind=which(allnm %in% nbdnm)
nb.rr=round(rr[,ind],digits=4)
colnames(nb.rr)=meths
rownames(nb.rr)=scens
nb.rr

grid.arrange(main="AERR",
textGrob("Poisson", y=0.2),
tableGrob(pois.rr),
textGrob("QuasiPoisson", y=0.2),
tableGrob(qpois.rr),
textGrob("Negative Binomial", y=0.2),
tableGrob(nb.rr), ncol=1
)


#active rates
ar.j2r=rbind(mean(poster1$mr.j2r.t),
	       mean(poster2$mr.j2r.t),
	       mean(poster3$mr.j2r.t),
	       mean(poster4$mr.j2r.t),
	       mean(poster5$mr.j2r.t))
ar.del=rbind(mean(poster1$mr.delta.t),
	   mean(poster2$mr.delta.t),
	   mean(poster3$mr.delta.t),
	   mean(poster4$mr.delta.t),
	   mean(poster5$mr.delta.t))

#placebo rates
pr.j2r=rbind(mean(poster1$mr.j2r.p),
	       mean(poster2$mr.j2r.p),
	       mean(poster3$mr.j2r.p),
	       mean(poster4$mr.j2r.p),
	       mean(poster5$mr.j2r.p))
pr.del=rbind(mean(poster1$mr.delta.p),
	   mean(poster2$mr.delta.p),
	   mean(poster3$mr.delta.p),
	   mean(poster4$mr.delta.p),
	   mean(poster5$mr.delta.p))
allnm=names(ar.j2r)

nbdnm=c("co_powernb", "dl_powernb", "j2r_powernb", "del_powernb" )
ind=which(allnm %in% nbdnm)
nb.rr=round(rr[,ind],digits=4)
colnames(nb.rr)=meths
rownames(nb.rr)=scens
nb.rr




nbdaerr1=data.frame(AERR=mar2_type1$rratio$dl_nbd.rratio,DOrate="16% NBD")
nbdaerr2=data.frame(AERR=mar2_type1b$rratio$dl_nbd.rratio,DOrate="34% NBD")
nbdaerr3=data.frame(AERR=mar2_type1$rratio$dl_quasipoisson.rratio,DOrate="16% QP")
nbdaerr4=data.frame(AERR=mar2_type1b$rratio$dl_quasipoisson.rratio,DOrate="34% QP")

nbdaerr=rbind(nbdaerr1,nbdaerr2,nbdaerr3,nbdaerr4)
m <- ggplot(nbdaerr, aes(x=AERR))
m +  geom_histogram(aes(y = ..density..),colour="grey80",fill="grey50") + geom_density() + facet_wrap(~DOrate)


nbdpval1=data.frame(pval=mar2_type1$pval$dl_nbd.pval,DOrate="16% NBD")
nbdpval2=data.frame(pval=mar2_type1b$pval$dl_nbd.pval,DOrate="34% NBD")
nbdpval3=data.frame(pval=mar2_type1$pval$dl_quasipoisson.pval,DOrate="16% QP")
nbdpval4=data.frame(pval=mar2_type1b$pval$dl_quasipoisson.pval,DOrate="34% QP")
nbdpval=rbind(nbdpval1,nbdpval2,nbdpval3,nbdpval4)
m <- ggplot(nbdpval, aes(x=pval))
m +  geom_histogram(aes(y = ..density..),colour="grey80",fill="grey50") + geom_density() + facet_wrap(~DOrate)


