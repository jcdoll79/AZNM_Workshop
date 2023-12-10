##########################################
#######Day 2.5: Introduction to R#########
##########################################


####Topics################################
#   R Studio and working with the console
#   Comments
#   Operators
#   Data Types
#   Creating New Variables
#   Complex data types
#   Functions
#   Packages
#   Coding Tips
#   Getting Help
#   Data exploration
#   Summarizing data
##########################################


#Comments#################################

#Just how do you add comments? 
#Comments are made by using the hash tag, #, character in front of the text. 
#The hash tag ONLY comments the single line.

#This is a comment line and will not be executed by R.
#Comments are useful to describe what you are trying to do



#Operators###############################

#Arithmetic operators
# + Addition
# - Subtraction
# * Multiplication
# / Division
# ^ or ** exponential

5+5
#COMMENT
45*2

#Logical operators
# > greater than
# >= greater than or equal to
# < less than
# <= less than or equal to
# == exactly equal to
# != Not equal to
5==5
5!=5
4!=5

#Other operators
# <- or = assignment operators



#Data Types#############################

#Numeric
#Assign the value of 20.5 to the object x
x <- 20.5
#print the values within the object x
x
#Determine the data type of object x
class(x)


#Integer
#Assign the integer value of 6 to the object y
y <- as.integer(6)
y2 <- 6
#print the values within the object y
y
#Determine the data type of object y
class(y)


#Using numeric variables 
#subtract y from x
x-y

#Character
#Assign the character string "Goodson Pond" 
#to the object z
z <- "Goodson Pond"
#print the values within the object z
z
#Determine the data type of object z
class(z)

#You can also convert numeric or 
#integer data to characters using a 
#special function "as.character()"
#Assign the character string "5.67" to 
#the object v
v <- as.character(5.67)
#print the values within the object v
v
#Determine the data type of object v
class(v)



#Creating new Variables/Objects########
#create a new object called "temp" and 
#assign it the value of 18.2
temp <- 18.2
#print the value stored in the object "temp"
temp



#Complex data types########################

#Vector
#A vector is a sequence of data of the *SAME* basic data type, 
#data types can't be mixed. A single vector can contain all numeric, 
#all integer, or all character types.

#create a vector of largemouth bass catch rates by lake. Each number represents 
#the total number of largemouth bass collected at a single lake. The name used 
#for the object which contains the total catch is "lmb".
lmb <- c(12,35,66,10,11,8)

#accessing elements within the lmb vector by specifying the index within 
#brackets. The following line returns the catch rate at the third lake.
lmb[5]

#create a vector of character data
char1 <- c("one","two","four","two","five","two")
char1 <- c("one","two","four","two","five","two")
char1[2]

#Matrix
#Two-dimensional layout of data. All data 
#must be the same type in a matrix.

lmb_byrow <- matrix(c(12,16,25,33,
                      35,40,10,45,
                      66,85,33,21,
                      10,10,25,10,
                      11,3,6,18,
                      8,7,12,22),  #the data elements entered one row at a time
                    nrow = 6,      #number of rows
                    ncol = 4,      #number of columns
                    byrow = TRUE)  #fill in the matrix by ROW

lmb_byrow #print the matrix in the console

#Return the dimensions of the lmb_byrow matrix
dim(lmb_byrow)

#Elements of a matrix are accessed by specifying 
#the row then column in brackets. To access 
#the catch from the third lake in the fourth 
#year use:
lmb_byrow[3,4]
#Access all catches from the third lake (row) 
#in the matrix, leave the column index empty
lmb_byrow[3,]
#Access all catches from the fourth year (column) 
#in the matrix, leave the row index empty
lmb_byrow[,3]

#This is the same data above but the data are 
#entered by column to demonstrate a different 
#way to entering data
lmb_bycolumn <- matrix(c(12,35,66,10,11,8,
                         16,40,85,10,3,7,
                         25,10,33,25,6,12,
                         33,45,21,10,18,22), #the data elements,one column at a time
                       nrow = 6,      #number of rows
                       ncol = 4,      #number of columns
                       byrow = FALSE)  #fill in the matrix by COLUMN
lmb_bycolumn #print the matrix in the console



#Data frame
#Another two-dimensional object but can contain different data types in 
#separate columns. For example, one column can be all numeric and 
#another column all characters

#Create a data frame with a column for lake names 
#and four years of data.
lmb_df <- data.frame("lake" = c("A","B","C","D","E","F"),
                     "2018" = c(12,35,66,10,11,8),
                     "2019" = c(16,40,85,10,3,7),
                     "2020" = c(25,10,33,25,6,12),
                     "2021" = c(33,45,21,10,18,22))

#NOTE column headers can't start with a number
#R will automatically add an "X" to "2018"

lmb_df

#Elements of a data frame can be accessed the same way as matrices.
#Catch from lake C in 2018
lmb_df[3,2]
#All catches from year 2018
lmb_df[,2]
#You can also return all elements from a single column by specifying the column name
lmb_df$X2018
lmb_df$X2018

lmb_df$X2018[2:4]

#List

#A list is a vector containing other objects. A list can contain multiple vectors, 
#multiple matrices, multiple data frames, or any combination. 
#The example below combined one vector, one matrix, and one data frame together 
#in a single list.

#See explanation of code for vector, matrix, and data frame above.
lmb <- c(12,35,66,10,11,8)
lmb_byrow <- matrix(c(12,16,25,33,
                      35,40,10,45,
                      66,85,33,21,
                      10,10,25,10,
                      11,3,6,18,
                      8,7,12,22),
                    nrow = 6,
                    ncol = 4,
                    byrow = TRUE)
lmb_df <- data.frame("lake" = c("A","B","C","D","E","F"),
                     "2018" = c(12,35,66,10,11,8),
                     "2019" = c(16,40,85,10,3,7),
                     "2020" = c(25,10,33,25,6,12),
                     "2021" = c(33,45,21,10,18,22))

#Combine the three objects created above into a 
#single list
lmb_list <- list(lmb, lmb_byrow, lmb_df)

#Print all elements from the list
lmb_list

#Access a member of the list using double brackets [[]] returns the matrix 
#from the lmb_list
lmb_list[[2]]

#Try to return the vector and data frame from the 
#list on your own.



#Import data file
Monroe_11 <- read.csv("Monroe_11.csv")
Monroe_11 <- read.csv("~/RWorkshop/AZNM_AFS/Monroe_11.csv")

#Break?

#Functions#######################################################

#A function performs a calculation or action to the data. 
#Nothing happens in R without using a function. You have already been using 
#functions in R! 
#data.frame() is a function

#A function includes the function name followed by parentheses. 
#The example below will calculate the mean of a specific row in a dataframe

#First, lets create a data frame with some data
lmb_df <- data.frame("lake" = c("A","B","C","D","E","F"),
                     "2018" = c(12,35,66,10,11,8),
                     "2019" = c(16,40,85,10,3,7),
                     "2020" = c(25,10,33,25,6,12),
                     "2021" = c(33,45,21,10,18,22))

lmb_df

#Now calculate the average catch across all lakes in a single
#Determine the average catch in 2018
mean(lmb_df[,2])

#Or reference the column header
mean(lmb_df$X2018)

#Determine the average catch at lake B
#Note that you must use a different function for row mean than column means.
#Specify the row and columns since the first columns includes lake
rowMeans(lmb_df[2,2:5])

#You can also save the output of a function into a 
#new object
Mean_2018 <- mean(lmb_df[,2])
#Print the mean from 2018 in the console
Mean_2018



#Packages#########################################################

#Packages contain the functions you will use to manipulate and analyze data.
#R comes pre-installed with several important packages.You will eventually 
#need to install other packages for specific analysis.
#The FSA package that we will use tomorrow does not come with R

#You only have to install a package once but you will need to update it periodically. 

#You can either use the "Install" button on the "Packages" tab,
#or use the following code
install.packages("tidyr")

#After installing a package, and before you use it, you have to tell R to load it.
#load the tidyr package. You will have to do this every time you open R.
library(tidyr)


#Getting Help###################################################

help(c)

#Or Google

#Exercise 1

#Data exploration/summary in R#########################################

#You might have to install the tidyr package
library(tidyr)

#Start with a data set but this time make it a tibble instead of a data frame
#A tibble won't add the annoying "X" to the numeric column header
lmb_df <- tibble("lake" = c("A","B","C","D","E","F"),
                 "2018" = c(12,35,66,10,11,8),
                 "2019" = c(16,40,85,10,3,7),
                 "2020" = c(25,10,33,25,6,12),
                 "2021" = c(33,45,21,10,18,22))


#Wide vs long format
#The lmb_df data frame is in in the "wide" format. 
#Wide format does not have values that repeat in any column. 
#Most R functions prefer data to be in the long format.

#This might be a change for users of other stats programs, such as SPSS. 
#SPSS likes data to be in the wide format. Long format contains values that repeat 
#in a column, usually the first column.

#There are many different ways to convert from wide to long or long to wide
#We will use the tidyr package

#pivot_longer() converts from wide to long
#pivot_wider() converts from long to wide

#Convert from wide to long
#Preview the data to remind us what we are looking at
lmb_df

#We will first introduce the "pipe" %>% symbol. The %>% is a useful operator to 
#chain multiple functions together. The basic syntax is: 
#data %>% function

#use Ctrl + Shift + M to add a pipe

#The object on the left of %>% becomes the first argument of the function to the right. The first argument of a 
#function is typically the data object that you want to apply the function to. Chained %>% can be used so the 
#output of one function is the input data for a new function:

#  -------> -------> -------> -------> ------->
#data %>% function1 %>% function2 %>% function3

#pivot_longer(cols = columns to pivot into longer,
#             names_to = the  name of the new column with groups,
#             values_to = the name of the new column that holds the values)


lmb_df %>% pivot_longer(cols = c("2018", "2019", "2020", "2021"),
                        names_to = "year",
                        values_to = "cpue")


#Now create a new data frame with the "long" data format
lmb_long <- lmb_df %>% pivot_longer(cols = c("2018", "2019", "2020", "2021"),
                                    names_to = "year",
                                    values_to = "cpue")


lmb_long

#Next we will convert our long format back to wide
#pivot_wider(names_from = column with the names of the groups,
#            values_from = column with the values for each group)


lmb_wide <- lmb_long %>% pivot_wider(names_from=year,
                                     values_from=cpue)

#There are many more arguments for the pivot_longer and pivot_wider functions. These allow for more complex 
#conversions such as transforming data and omitting data. 
#We won't go into those details here but they do exist!


#The next set of functions will use a data set from the FSAdata package. You might have to install the 
#FSAdata package.

#Load the FSAdata package
library(FSAdata)

#Load the bass data set from Florida
BassFL <- FSAdata::BassFL
#Check out the help file associated with it to see 
#what the columns represent
help("BassFL")


#Exploring the BassFL data

#There are times you need to know the list of values in a column. This can be helpful when you need to filter 
#data (which we'll do later).

#What are the unique species?
unique(BassFL$age)
#what are the unique locations? Try on your own.


#What are the minimum and maximum number caught across species, location, and year?
range(BassFL$num)
#What about the minimum and maximum age? Try on your own
range(BassFL$age)


#Next, we will demonstrate filter, sorting, summarizing, 
#and adding new columns

names(BassFL)
head(BassFL)
unique(BassFL$species)
#Filter all records of Suwanee bass
#Note there are a few different packages that use the "filter" function. We specifically want the "filter" 
#function from the "dplyr" package. To do this we need to specify 
#the package then :: then function 
#package::function
Suwanee <- BassFL %>% 
           dplyr::filter(species == "Suwanee")
Suwanee


#Filter only Suwanee bass that are age 4 and older
Suwanee_4p <- BassFL %>%
              dplyr::filter(species == "Suwanee" & age >= 4)
Suwanee_4p

#Using the code above as a template, modify to filter 
#for all observations from SantaFe that captured more 
#than 39 fish

#Answer###############
Suwanee_n40 <- BassFL %>%
               dplyr::filter(loc == "SantaFe" & num >= 40)
Suwanee_n40
#End Answer##################

#Sort##########
#The arrange function sorts alphabetically and ascending by default
Bass_sort_spnum <- BassFL %>%
                   dplyr::arrange(species,num)
Bass_sort_spnum

#You can specify descending using desc() function
Bass_sort_spnum_desc <- BassFL %>%
                        dplyr::arrange(desc(species),desc(num))
  

#Summarizing data#################

#Determine the average age by species and location
#To accomplish this, we will chain together multiple functions using pipes.

BassFL %>%    #specify data set
  dplyr::group_by(species,loc) %>% #identify groups to summarize
  dplyr::summarise(Mean_age = mean (age, na.rm = TRUE))
  
#Use what you have learned above to
#1) filter for only Suwanee bass
#2) group data by location
#3) create a new column titled "Mean_num" that will include the average catch by location
#4) Sort the data by ascending "Mean_num"


#Answer###########
BassFL %>%
  dplyr::filter(species == "Suwanee") %>%
  dplyr::group_by(loc) %>%
  dplyr::summarise(Mean_num = mean (num, na.rm = TRUE)) %>%
  dplyr::arrange(Mean_num)

#End Answer#######


#Count the number of observations per species#########
BassFL %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(num_obs = dplyr::n())



#Exercise 2

#Basic graphing in R#########################
#Base R graphs are helpful to quickly visualize data but they are not very aesthetically pleasing.
#We will use ggplot2 instead. Advanced graphs and ultimate flexibility can be achieved using advanced 
#ggplot2 functions. We won't cover advanced ggplot2 here - that could be an entire day!

#You might have to install the ggplot2 package. 
#Install ggplot2 using the "Install" button under the 
#"Packages" tab. 

#load ggplot2
library(ggplot2)
library(dplyr)
library(FSAdata)

#Scatterplot#########

#Create a scatterplot of number at age
#Filter the bass data to only include Largemouth from SantaFe in 2001
LMB_SantaFe <- BassFL %>% 
  dplyr::filter(species == "Largemouth" & 
                loc == "SantaFe" & 
                year == 2001)

#ggplot2 uses "layers" to build a figure. All ggplot2 
#figures begin with specifying the data

#ggplot(data, aes(x = x values, y = y values))

#Produces a blank graph
ggplot( LMB_SantaFe, aes(x = age, y = num) ) 

#After declaring the data layer, you can add elements to build it using the "+" symbol
#The geom_point function tells R to make a scatterplot of points using the specified data
ggplot(LMB_SantaFe, aes(x = age, y = num)) +
  geom_point() 

#Change the point size and color
ggplot(LMB_SantaFe, aes(x = age, y = num)) +
  geom_point(shape = 1, size = 10, color = "blue") 
#Try different values for size, shape, and different 
#colors.



#Histogram#######
#Load the Bonito data set from the FSAdata package
Bonito <- FSAdata::Bonito
help(Bonito)

#Allow ggplot to pick bin width
ggplot(Bonito, aes(x = fl)) + 
  geom_histogram()

#Adjust bin width
ggplot(Bonito, aes(x=fl)) + 
  geom_histogram(bins = 20)

#Adjust the color of the bars
#fill will change the fill color, color will change 
#the bar borders
ggplot(Bonito, aes(x=fl)) + 
  geom_histogram(bins = 20, fill = "pink", color = "black")


#Boxplot#############
#Boxplot of fork length by sex
ggplot(Bonito, aes(x = sex, y = fl)) +
  geom_boxplot()


#Overlay raw data that are "jittered"
ggplot(Bonito, aes(x = sex, y = fl)) +
  geom_boxplot() +
  geom_jitter()


#Modify axes#########
#Controlling the axes format of the scatterplot introduced earlier
#Modify axes labels
ggplot(LMB_SantaFe, aes(x = age, y = num)) +
  geom_point(shape = 8, size = 4, color = "darkgreen") +
  xlab("Age (years)") + 
  ylab("Number")


#Change the font size
ggplot(LMB_SantaFe, aes(x = age, y = num)) +
  geom_point(shape = 8, size = 4, color = "darkgreen") +
  xlab("Age (years)") + 
  ylab("Number") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

