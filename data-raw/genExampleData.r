#creates simulated dataset used in first part of vignette
set.seed(698123)

#sample size
n <- 500
#generate randomised treatment group
z <- 1*(runif(n)<0.5)

#specify true log rate ratio for effect of treatment
trtEffect <- -1

#simulate subject level random effect
k <- 4
gamma <- rgamma(n,shape=k,scale=1/k)
rate <- gamma*exp(trtEffect*z)

#planned follow-up is 1 year. Simulate random time around 1 year for each
#subject, which would be their follow-up if they completed the study.
fullFupTime <- 1+0.1*runif(n)

#simulate time at which they might drop out, between 0 and 1
dropOutTime <- runif(n)

#simulate number of events up to dropOutTime
y1 <- rpois(n, lambda=rate*dropOutTime)
#simulate number of events in remaining time, assuming patient remains
#on their randomised treatment
y2 <- rpois(n, lambda=rate*(fullFupTime-dropOutTime))
#generate number of patients that would be seen if patient was on placebo
#in the period from dropOutTime to fullFupTime
y2placebo <- rpois(n, lambda=gamma*exp(0)*(fullFupTime-dropOutTime))

#model dropout process, with dropout depending on count up to dropOutTime
do_xb <- -1+1*log(y1+1)
do_pr <- exp(do_xb)/(1+exp(do_xb))
do <- 1*(runif(n)<do_pr)

#we will assume patients who dropout in active arm switch to placebo,
#so set their count in the second period to that simulated under placebo
y2[(do==1) & (z==1)] <- y2placebo[(do==1) & (z==1)]

#generate observed count for each patient
y <- y1+y2*(1-do)

#generate a variable containing each patient's actual follow-up time
fupTime <- fullFupTime
fupTime[do==1] <- dropOutTime[do==1]

#create data frame to save with package
simData <- data.frame(z,y,fupTime)

devtools::use_data(simData, overwrite=TRUE)