#Exercise 3 Answer Key
library(FSA)
library(ggplot2)
library(tidyr)
library(dplyr)

#1.	Generate a vector of percentage of fish in 10mm length intervals using the 
#“BluegillLM” data set from the FSAdata package.
BluegillLM <- FSAdata::BluegillLM 

#Add a length category column
BluegillLM$lcat10 <- lencat(BluegillLM$tl,w=10)

#Check data to confirm
headtail(BluegillLM)

#Create a length frequency table
#Using the xtabs function, requires the column name 
#with length category and the dataset
BLGFreq10 <- xtabs(~lcat10, data = BluegillLM)
#print frequency table
BLGFreq10

#Percentage of fish in each interval using the 
#prop.table() function
BLGPer10 <- round( prop.table(BLGFreq10) * 100 , 1)
#print proportion table
BLGPer10


#2.	Create a length frequency histogram of total length (tl column) using the “BluegillLM” 
#data set from the FSAdata package.

#A histogram can be useful to visualize 
#length frequencies
ggplot(BluegillLM, aes(x = tl)) + 
  geom_histogram(breaks = seq(from = 50 , to = 230, by = 10),
                 fill = "orange", color = "black")

#3.	Create a frequency table of Bluegill PSD size groups using the “BluegillLM” 
#data set from the FSAdata package. 


#Assign length category to each fish to the 
#Largemouth Bass data frame
#Load the BluegillLM data set
BluegillLM  <- FSAdata::BluegillLM 

#First we need to pull out the Largemouth Bass 
#length categories
blg.cuts <- psdVal("BLuegill")

#Filter based on tl greater then stock length.
#This eliminates fish less than stock length
#and add the length category column
BLG_SS <- BluegillLM  %>%
  filter( tl >= blg.cuts["stock"]) %>%
  mutate( gcat = lencat(tl, breaks = blg.cuts,
                        use.names = TRUE))

headtail(BLG_SS)

#Calculate a frequency table across size groups
xtabs(~gcat, data = BLG_SS)



#4.	Calculate the PSD-Q and PSD-P of Bluegill using the BluegillLM dataset 
#from the FSAdata package.

#Calculate all PSD-X values and 95% Confidence 
#Intervals. Requires length column, data set, 
#species, and what type. Types available are 
#incremental (e.g., Stock to Quality) or
#traditional (e.g., Quality, Preferred, Memorable)
psdCalc(~tl, data=BLG_SS, 
        species = "Bluegill", 
        what = "traditional")


#5.	Create a scatterplot of Bluegill total length (x-axis) and weight (y-axis) on the 
#natural log scale using the “BluegillLM” data set from the FSAdata package.

#Weight-length relationships#######
#Load Chinook data from the FSA package
BluegillLM <- FSAdata::BluegillLM

#Plot TL vs TW on the natural log scale
ggplot(BluegillLM, aes(x = log(tl), y = log(wght))) +
  geom_point() +
  xlab("ln(total Length) (mm)") + 
  ylab("ln(total Weight) (g)")


#6.	Estimate coefficients of a weight-length model using the “BluegillLM” 
#data set from the FSAdata package. Also determine the 95% confidence intervals for the 
#intercept and slope.

#Estimate parameters of the weight-length model 
#using lm()
lm1 <- lm(log(wght) ~ log(tl), data=BluegillLM)

#Extract summary information and send to a new object
sumlw <- summary(lm1)
sumlw

#Extract coefficients and send to a new object
coeflw <- coef(lm1)
coeflw

#Extract confidence intervals for coefficients
confinlw <- confint(lm1)
confinlw

