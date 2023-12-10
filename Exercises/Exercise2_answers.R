#Exercise 2 Answer Key
library(tidyr)
library(ggplot2)

#1.	Load the “Iris data” using the following code
data("iris")

#2.	Convert the wide iris data to long. 
#   a.	Columns to convert to long include: Sepal.Length, Sepal.Width, Petal.Length, and Petal.Width
#   b.	Name the new names column “measurement”
#   c.	Name the new value column “length”

iris_long <- iris %>% pivot_longer(cols = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
                                    names_to = "measurement",
                                    values_to = "length")

#3.	Return a vector of unique species names in the “iris_long” data frame

unique(iris_long$Species)

#4.	Using the “iris_long” data frame, filter only Sepal Lengths that are greater than 5.0 and 
#assign the results to a new data frame named “filtered_iris’
filtered_iris <- iris_long %>%
                dplyr::filter(measurement == "Sepal.Length" & length > 5.0)
filtered_iris

#5.	Using the "iris_long" data frame, select only the setosa species, group the data by measurement, 
#calculate the mean length for each measurement, and sort by the mean length for each measurement.

iris_long %>%
  dplyr::filter(Species == "setosa") %>%
  dplyr::group_by(measurement) %>%
  dplyr::summarise(Mean_length = mean (length, na.rm = TRUE)) %>%
  dplyr::arrange(Mean_length)


#6.	Create a scatterplot of setosa Sepal. Length (x-axis) and Sepal.Width (y-axis). 
#Note you will have to filter the original iris data set.

irissub <- iris %>% 
  dplyr::filter(Species == "setosa")

ggplot(irissub, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point() 


#7.	Adjust the scatterplot created in 6 to the following: use the shape number 9, 
#increase the size to 3, change the color to anything you would like, 
#change the x-axis label to "Sepal length (cm)" and the x-aix label to "Sepal width (cm).

ggplot(irissub, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(shape = 9, size = 3, color = "darkgreen") +
  xlab("Sepal length (cm)") + 
  ylab("Sepal width (cm)")
