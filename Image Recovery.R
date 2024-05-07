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
image(image, axes = FALSE)

# Add Gaussian noise to the image
set.seed(1234) #set a seed for reproducibility
noisy_image <- image + rnorm(length(image), mean = 0, sd = 1)  # Adjust mean and sd as needed

#Display the noisy spiral image. Note that the spiral is barely visible.
image(noisy_image, axes=FALSE)

######################################
#Process the noisy image to recover the original.
#I start by renaming image to X, for simplicity
X <- image

#Singular Value decomposition
svdX <- svd(X)
si <- seq(1, length(svdX$d), 1)

plot(si, svdX$d, main="Singular Values of Noisy Matrix", xlab="Index of Singular Value", 
     ylab="Singular Value")

energy <- svdX$d/sum(svdX$d)
cumenergy <- cumsum(energy)
plot(si, cumenergy, main="Cumulative Energy of Singular Values", xlab="Index of Singular Value",
     ylab="Cumulative Energy of Singular Values")


#For comparison, use the full SVD for image reconstruction
full <- svdX$u[,1:300] %*% diag(svdX$d[1:300], 300) %*% t(svdX$v[,1:300])
image(full, axes=FALSE)

#Capture 90% Energy
ninety <- which(cumenergy>=0.9)[1]
u90 <- svdX$u[,1:ninety] 
v90 <- svdX$v[,1:ninety]
d90 <- svdX$d[1:ninety]
M90 <- diag(d90, ninety)
image90 <- u90 %*% M90 %*% t(v90)

#Now, plot the image with 90% energy recovery
image(image90, main = "90% Energy")

#calculate error
(error <- tr(t(image90-image)%*%(image90-image))/300^2) #error is based on frobenious norm.


#davish and donoho with sigma estimate
k <- 0.4389
w <- (4/sqrt(3))/k
svd4 <- which(svdX$d >= w)
usvd4 <- svdX$u[, 1:4]
vsvd4 <- svdX$v[, 1:4]
dsvd4 <- svdX$d[1:4]
Msvd4 <- diag(dsvd4, 4)
imagesvd4 <- usvd4 %*% Msvd4 %*% t(vsvd4)
image(imagesvd4, main="tSVD with unknown sigma")

#caluclate error
tr(t(imagesvd4-image)%*%(imagesvd4-image))/300^2

#calculate error
(error <- tr(t(image90-image)%*%(image90-image))/300^2)


#just for fun, see what 80% energy looks like.
eighty <- which(cumenergy>=0.8)[1]
u80 <- svdX$u[,1:eighty] 
v80 <- svdX$v[,1:eighty]
d80 <- svdX$d[1:eighty]
M80 <- diag(d90, eighty)
image80 <- u80 %*% M80 %*% t(v80)

#Now, plot the image with 90% energy recovery
image(image80, main = "80% Energy")
