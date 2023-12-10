#Exercise 1 Answer Key

#1.	Assign the value of 120 to x
x <- 120

#2.	Assign the value of 345 to y
y <- 345

#3.	Create a new variable z with the value y – x
z <- y-x

#4.	Return the value of z to the console
z

#5.	Create a matrix named “dat1” that has five row and four columns with the following data:
#      266	369	452	555
#      245	456	444	576
#      221	447	489	589
#      243	299	465	512
#      285	325	477	499

dat1 <- matrix(c(266,	369,	452,	555,
                 245,	456,	444,	576,
                 221,	447,	489,	589,
                 243,	299,	465,	512,
                 285,	325,	477,	499), #the data elements,one column at a time
               nrow = 5,      #number of rows
               ncol = 4,      #number of columns
               byrow = TRUE)  #fill in the matrix by ROW

#6.	 Determine the dimensions of the “dat1” matrix using a function.
dim(dat1)

#7.	Return the second row in the “dat1” matrix created in number 5
dat1[2,]

#8.	Create a data frame named “Lake1” that has six row and two columns with the following data:
#  Species	TL
#   LMB	195
#   LMB	210
#   LMB	222
#   LMB	168
#   BLG	95
#   BLG	125

Lake1 <- data.frame (Species = c("LMB", "LMB", "LMB", "LMB", "BLG", "BLG"),
                     TL = c(195, 210, 222, 168, 95, 125))


#9.	Return the second column (total length column) of the “Lake1” data from two different ways.
Lake1[,2]      
Lake1$TL

#10. Install the Hmisc package.