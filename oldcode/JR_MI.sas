
%let nimpute=10;

libname library "\\emea.astrazeneca.net\sweden\Molndal\Users 05\kdss481\Documents\Training\Missing data\INR\data";

proc import 
datafile="\\emea.astrazeneca.net\sweden\Molndal\Users 05\kdss481\Documents\Training\Missing data\INR\data\mardat1.csv"
out=library.simN replace;
run;

data DS0;
	set library.simN;
* Log study time for use as offset;
	offset1=log(Study_Time);
* Add record number (this is used for merging data back together;
	Record=id;
* Maxtime should be same throughout and is the time which imputed number of events is up to (i.e. length of study;
	maxtime=1;
run;

* Function to do the calculations in the paper;
proc fcmp outlib=work.funcs.trial;
function extra(count,dispersion, predict1, predict2, time, maxtime);
	invk=1/dispersion;
	b1t1=exp(predict1)*Time;                     /*mean number of exacerbations before withdrawal*/
	b2t2=exp(predict2)*max((MaxTime-Time),0);     /*mean number of exacerbations after withdrawal*/
	p2=(invk+b1t1)/(invk+b1t1+b2t2);
* Event is not possible or so unlikley that we can ignore;
	if  p2 > 0.999999 then return(Count);
* InvK is so large that we use Poisson (Dispersion is very cloe to zero);
	else IF invK > 1.E20 then return(Count+ rand('Poisson',b2t2));
	else return(Count+ rand('NEGATIVEBINOMIAL',p2, invK+Count));
endsub;

function getrate(count,dispersion, predict1, predict2, time, maxtime);
	invk=1/dispersion;
	b1t1=exp(predict1)*Time;                     /*mean number of exacerbations before withdrawal*/
	b2t2=exp(predict2)*max((MaxTime-Time),0);     /*mean number of exacerbations after withdrawal*/
	p2=(invk+b1t1)/(invk+b1t1+b2t2);
	return(p2*(invK+Count));
	/*Rate=k*rate*((invk + obs))/(1+k*rate)) */
endsub;
quit; 

data start;run;

data start;run;
options cmplib=work.funcs;

%macro simMI(inset,N);
%let i = 1;

%DO i = 1 %TO &N;
data seti;
	set &inset (where=(sim eq &i));
run;

* Build sample from Bayesian posterior for parameters;
ODS RESULTS OFF;
ODS select none;
proc genmod data=seti;
	class trt;
	model count = trt /noint link=log dist=negbin offset=offset1 ;
	bayes OUTPOST=OUTP seed=12345 thin=20 nmc=%eval(&nimpute*20) NBI=1000; 
run;
ods Results on;
ODS select all;

* Join the two data sets;
proc sql;
create table jr as
select
	A.*, iteration, Dispersion, TrtActive, trtPlacebo
	from seti A , outp ;
quit;

proc sort data= jr; by iteration record;
* Here is the imputation stage;
* Note that INVK is the value of k in the paper;
data MI;
	set jr;
	by iteration record;
/*Seed for the RAND function*/
	call streaminit(457);
    mup=exp(trtPlacebo);
	mua=exp(trtActive);
* MAR;
	method="MAR";
	*if trt= 0 then predict1=trtPlacebo;
	*if trt= 1 then predict1=trtTreatment;
	if trt= "Placebo" then predict1=trtPlacebo;
	if trt= "Active" then predict1=trtActive;
	predict2=predict1;
	imputed=extra(count,dispersion,predict1,predict2,Study_time,maxtime);
	postwdr=getrate(count,dispersion,predict1,predict2,Study_time,maxtime);     /*rate used in imputations*/
	output;

* J2R;
	method="J2R";	
	predict2=TrtPlacebo;
	imputed=extra(count,dispersion,predict1,predict2,Study_time,maxtime);
	postwdr=getrate(count,dispersion,predict1,predict2,Study_time,maxtime);     /*rate used in imputations*/
	output;

* CR;
	method="CR";
	predict1=TrtPlacebo;
	imputed=extra(count,dispersion,predict1,predict2,Study_time,maxtime);
	postwdr=getrate(count,dispersion,predict1,predict2,Study_time,maxtime);     /*rate used in imputations*/
	output;
run;

proc sort data=mi;
	by method iteration trt;
run;

proc means data=mi noprint;
	by method iteration trt;
	var count cc imputed study_time;
	output out=events sum=events ev_c ev_imp time;
run;

data events;
	set events;
	ev_exp=round(events*(_freq_/time));
	sim=&i;
run;

data eventsN; set eventsN events; run;

%put &i;
%end;
%mend;

data eventsN;run;
%simMI(DS0,5);

data eventsN; set eventsN; 
	rate= events/time;
	rate_imp= ev_imp/_freq_;
	rate_c= ev_c/_freq_;
run;

proc sort data=eventsN; by method sim trt;
proc means data=eventsn noprint;
	where sim ne .;
	by method sim trt;
	var events ev_imp ev_exp ev_c rate rate_imp rate_c;
	output out = evmeans mean=events ev_imp ev_exp ev_c rate rate_imp rate_c;
run;

proc sort data=evmeans; by method sim trt rate_c;
proc transpose data=evmeans out=evmeanst;
	by method sim trt rate_c;
	var rate_imp;
run;

data mi;
set mi;
 psi=0.8/1;
 ptst=psi/(1+psi);
 tst=rand('NEGATIVEBINOMIAL', 1-ptst, 1);
run;
 proc means data=mi;
 	var tst;
 run;

proc export 
outfile="\\emea.astrazeneca.net\sweden\Molndal\Users 05\kdss481\Documents\Training\Missing data\INR\data\sasres.csv"
data=evmeanst replace;
run;


proc import 
datafile="\\emea.astrazeneca.net\sweden\Molndal\Users 05\kdss481\Documents\Training\Missing data\INR\data\mardat2.csv"
out=library.simN replace;
run;

data DS0;
	set library.simN;
* Log study time for use as offset;
	offset1=log(Study_Time);
* Add record number (this is used for merging data back together;
	Record=id;
* Maxtime should be same throughout and is the time which imputed number of events is up to (i.e. length of study;
	maxtime=1;
run;


data start;run;
data eventsN;run;
%simMI(DS0,100);

data eventsN; set eventsN; 
	rate= events/time;
	rate_imp= ev_imp/_freq_;
	rate_c= ev_c/_freq_;
run;

proc sort data=eventsN; by method sim trt;
proc means data=eventsn noprint;
	where sim ne .;
	by method sim trt;
	var events ev_imp ev_exp ev_c rate rate_imp rate_c;
	output out = evmeans mean=events ev_imp ev_exp ev_c rate rate_imp rate_c;
run;

proc sort data=evmeans; by method sim trt rate_c;
proc transpose data=evmeans out=evmeanst;
	by method sim trt rate_c;
	var rate_imp;
run;

proc export 
outfile="\\emea.astrazeneca.net\sweden\Molndal\Users 05\kdss481\Documents\Training\Missing data\INR\data\sasres2.csv"
data=evmeanst replace;
run;






ODS RESULTS OFF;
ODS select none;
option nonotes;
proc genmod data =MI;
class trt ;
model Imputed = trt  /Dist=Negbin link=log ;
lsmeans trt / diff ;
by method iteration;
ods output  diffs=LSMD ;
run;
option notes;
ods Results on;
ODS select all;


proc sort data=lsmd;
by method trt _trt;
run;

ODS RESULTS OFF;
ODS select none;
	proc mianalyze data=lsmd alpha=0.05;
	modeleffects estimate;
	STDERR Stderr;
	by method trt _trt;
	ods output parameterestimates=RES_D2;
	run;
ods Results on;
ODS select all;


data dsout;
set res_d2;
keep  Method Nimpute trt _trt estimate  invest probt  LowCL UppCL iLowCL iUppCL MCSE;
attrib Estimate label="Estimate" format=d12.3;
Nimpute=&Nimpute;
estimate=exp(estimate);
invest=1/estimate;
LowCL=exp(LCLMean);
iLowCL=1/Lowcl;
UppCL=exp(UCLMean);
iUppCL=1/UppCL;
output; 
run;

proc print data=dsout;
run;
