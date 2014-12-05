path="\\\\emea.astrazeneca.net\\sweden\\Molndal\\Users 05\\KDSS481\\Documents\\Training\\Missing data\\INR"

#Functions to generate data
#In the example, nb_complete, nb_incomplete and dropout_time4 functions are used
source(paste(path,"/Code/Latest/Robertscodes.R", sep=""))

#Functions to do multiple imputations
#In the example cim function is used
source(paste(path,"/Code/Latest/cim_v1.1.R", sep=""))

#Functions to evaluate missing data methods
#In the example sim_mar2 function is used
source(paste(path,"/Code/Latest/power_functionsv1.2.R", sep=""))

NS=1000
sasname=paste(path,"/Data/SASvsRparam3.csv",sep="")
# MAR (drop out rates with variation - under null, equal DO rates) - Type I Error 
mar2_type1np <- sim_mar2(n=228,rate.placebo=0.88,rr=0.4,T=1,
                      dispersion.placebo=0.9, dispersion.treatment=0.9,
                      rates=seq(0.001,40,2),drop_rate_variance=0.05,
                      N.rep=NS, N.mi=10, alpha=0.04, sasname)