#Exercise 4 Answer Key
library(FSA)
library(dplyr)
library(ggplot2)
library(nlstools)

#1.	Calculate standard weight and relative weight of Cisco using the “CiscoTL” 
#data set in the FSAdata package. Add these two columns to the CiscoTL data frame.
#Load the “CiscoTL” data set using:
CiscoTL <- FSAdata::CiscoTL

#Return a simplified object for calculation
wsCIS <- wsVal("Cisco", units=c("metric"), 
               simplify = TRUE)

#Add Ws and Wr column
CiscoTL <- CiscoTL %>%
  mutate(logW = log10(weight),
         logL = log10(length),
         Ws = 10 ^(wsCIS[["int"]] + wsCIS[["slope"]] * logL ),
         Wr = weight/Ws * 100)

#Note there are missing weight which will result in NA being
#returned for relative weights.


#2.	Use the following data to create a scatterplot of count (y-axis) vs 
#age (x-axis) and determine instantaneous total mortality (Z)

#Create a data frame for say, Brook Trout
bkt <- data.frame(age=1:8,
                  ct=c(74,210,165,92,82,50,25,10))
bkt

#Create a scatterplot with log(ct) to 
#identify the descending limb of the catch curve
ggplot(bkt, aes(x = age, y = log(ct))) +
  geom_point() +
  xlab("Age (years)") + 
  ylab("ln(count)")

bktcc <- catchCurve(ct ~ age, data = bkt, ages2use=2:8)

#The summary function will return
#the instantaneous mortality (Z)
#and annual mortality (A)
summary(bktcc)

#the plot() function will create a plot of the 
#catch curve
plot(bktcc)

#3.	Create a scatterplot and determine parameters of the von Bertalanffy growth model 
#with 95% confidence intervals using “Bonito” data set in the FSAdata package. 
#Bonito length is recorded as “fl”.

Bonito <- FSAdata::Bonito

ggplot(Bonito, aes(x = age, y = fl)) +
  geom_point() +
  xlab("Age (years)") + 
  ylab("Total Length (mm)")

vbT <- vbFuns("typical", simple=FALSE)
fitBonito <- nls(fl ~ vbT(age, Linf ,K, t0),
                 data=Bonito,
                 start=vbStarts(fl ~ age, 
                                data = Bonito, 
                                type="typical"))

sumBonito <- summary(fitBonito, correlation = TRUE)
sumBonito
coefBonito <- coef(fitBonito)
coefBonito
confinT <- confint2(fitBonito)
confinT


#4. Estimate parameters of the Ricker stock-recruitment model using the HerringBWE 
#data set from the FSAdata package
HerringSR <- FSAdata::HerringBWE

headtail(HerringSR)

#Plot spawners (stock) vs recruits
ggplot(HerringSR, aes(x = ssb, y = recruits)) +
  geom_point() 


#Ricker function
#E[R|S] = alpha * S * exp(-beta * S)

#We will use the nls function to fit this 
#non-linear model
#Requires starting values
svR <- srStarts(recruits ~ ssb, 
                data = HerringSR, 
                type = "Ricker")
svR

#Obtain Ricker function from FSA
rckr <- srFuns("Ricker")

#Fit Ricker function to stock and recruitment data
srR <- nls(recruits ~ rckr(ssb,a,b), 
           data = HerringSR, 
           start=svR)
summary(srR)

#Coefficients with 95% Confidence Intervals
cbind(estimates=coef(srR), confint(srR))

#5.	Create a scatterplot with predictions of the Ricker model fit in 4 above. 
#Note, select an appropriate break in the range of spawning stock biomass for predictions.
#Visualize the model fit
#Range of spawning stock
x <- seq(from=min(HerringSR$ssb), 
         to = max(HerringSR$ssb), 
         by= 10)

#Predict recruitment from model fit above
pR<- rckr(x, a=coef(srR))
#combine in a data frame
CombSR <- data.frame(x = x, pR = pR)

#Plot predictions with raw data
ggplot() +
  geom_point(aes(x = CombSR$x, y = CombSR$pR), 
             color = "orange", size = 3) +
  geom_point(aes(x = HerringSR$ssb, 
                 y = HerringSR$recruits)) +
  xlab("Spawning stock biomass") + 
  ylab("Predicted recruitment")
