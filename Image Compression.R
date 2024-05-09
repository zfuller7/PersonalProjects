#Image resolution improvement
#load libraries 
library(psych) #needed for tr() function

# Generate an image to start with.
width <- 300  # Adjust width and height as needed
height <- 300
center <- c(width / 2, height / 2)
theta <- seq(0, 6 * pi, length.out = 1000)  # Adjust the number of points for smoother or coarser spiral
radius <- seq(0, min(width, height) / 4, length.out = length(theta))

# Compute spiral coordinates
spiral_x <- center[1] + radius * sin(theta)
spiral_y <- center[2] + radius * cos(theta)

# Create an empty image
image <- matrix(0, nrow = height, ncol = width)

# Plot the spiral on the image
for (i in 1:length(spiral_x)) {
  x <- round(spiral_x[i])
  y <- round(spiral_y[i])
  if (x >= 1 && x <= width && y >= 1 && y <= height) {
    image[y, x] <- 1  # Set pixel value to 1 (white) for the spiral
  }
}

# Display the original spiral image
image(image, axes = FALSE, main="Original")



######################################
#Functions for image recovery
#This first function applies SVD to the image
image_svd <- function(input_image){
  image_decomposition <- svd(input_image)
  energy <- image_decomposition$d/sum(image_decomposition$d)
  cummulative_energy <- cumsum(energy)
  return(list(u=image_decomposition$u, v=image_decomposition$v, d=image_decomposition$d, 
              cummulative_energy=cummulative_energy))
}

#This function takes as input the object created by image_svd and optional energy_level.
image_reconstruction <- function(image_svd, energy_level=0.90){
  energy_filter <- which(image_svd$cummulative_energy >= energy_level)[1]
  u <- image_svd$u[,1:energy_filter]
  v <- image_svd$v[,1:energy_filter]
  d <- image_svd$d[1:energy_filter]
  M <- diag(d, energy_filter)
  image_return <- u %*% M %*% t(v)
  energy_level=energy_level
  return(list(image_return=image_return, energy_level=energy_level))
}

#This should be used when you want to set a number of singular values instead of energy level.
image_reconstruction2 <- function(image_svd, sv=100){
  u <- image_svd$u[,1:sv]
  v <- image_svd$v[,1:sv]
  d <- image_svd$d[1:sv]
  M <- diag(d, sv)
  image_return <- u %*% M %*% t(v)
  sv=sv
  return(list(image_return=image_return, sv=sv))
}

image_plot <- function(img){
  image(img$image_return, main=paste0("Image with ", img$energy_level*100, "% Energy"), axes=FALSE)
}

#works with image_reconstruction2
image_plot2 <- function(img){
  image(img$image_return, main="Image with Custom Cut-off", axes=FALSE)
}
###############################################
#This section contains graphs for explaining the "energy" part of the process.

#Singular Value decomposition
X <- image
svdX <- svd(X)
si <- seq(1, length(svdX$d), 1)

plot(si, svdX$d, main="Singular Values of Noisy Matrix", xlab="Index of Singular Value", 
     ylab="Singular Value")

energy <- svdX$d/sum(svdX$d)
cumenergy <- cumsum(energy)
plot(si, cumenergy, main="Cumulative Energy of Singular Values", xlab="Index of Singular Value",
     ylab="Cumulative Energy of Singular Values")
abline(h=0.9)
abline(v=125)
############################################################
#set parameters to plot all 4 graphs together
par(mfrow=c(2,2))

#Graph original image
image(image, main="Original", axes=FALSE)
#using visual cut-off from Cumulative Energy Graph
original_svd <- image_svd(image)
new_image <- image_reconstruction2(original_svd, sv=125)
image_plot2(new_image)

#calculate error
(error <- tr(t(new_image$image_return-image)%*%(new_image$image_return-image))/300^2) #error is based on frobenious norm.

#Capture 90% Energy
new_image <- image_reconstruction(original_svd, energy_level=0.9)
image_plot(new_image)

#calculate error
(error <- tr(t(new_image$image_return-image)%*%(new_image$image_return-image))/300^2) #error is based on frobenious norm.

#Capture 80% Energy
new_image <- image_reconstruction(original_svd, energy_level=0.8)
image_plot(new_image)

#calculate error
(error <- tr(t(new_img2$image_return-image)%*%(new_img2$image_return-image))/300^2)




