setwd("C:/Users/longle/")
rm(list=ls(all=TRUE))
# read data
# separation ","
x <- read.table("C:\\Users\\longle\\Dropbox\\data\\train3.data",sep=",")
x <- as.matrix(x)
dim(x)
# [1] 658 256


# plot the first image of 3
windows()
image(matrix(x[1,],16))
# gray scale
image(matrix(x[1,],16),col=gray((32:0)/32))
# upside down?

#Long's edit
#One picture: x2 = matrix(x[1,],16), rx2 <- x2[, 16:1]
#All pictures:
#rx <- x[, 241: 256]
#for (i in 15:1) {
#    rx <- cbind(rx, x[, (16*i + 1): (16*i + 16)])
#}

# convert to a 3D array
y <- array(x,c(nrow(x),16,16))
# correct the upside down
y <- y[,,16:1]
windows()
image(y[1,,],col=gray((32:0)/32))
# looks all right
# change the margin!!!
windows()
par(mfrow = c(4,4),mar=rep(1,4))
for(i in 1:16){
  image(y[i,,],col=gray((32:0)/32))
}
# convert back to 2D matrix
x <- matrix(y,ncol=256)


# compute the mean image
x.mean <- apply(x,2,mean)
# subtract the mean image !!!!!!!!!!
x.demean <- x-rep(1,nrow(x))%*%t(x.mean)
# SVD of the centered matrix
x.svd <- svd(x.demean)
# plot mean, the first and second eigen-image
windows(height=3,width=8)
par(mfrow = c(1,3))
image(matrix(x.mean,16),col=gray((32:0)/32))
image(matrix(-x.svd$v[,1],16),col=gray((32:0)/32))
image(matrix(x.svd$v[,2],16),col=gray((32:0)/32))

###########################################################################
# analyze 4 and 9 together
# load 4
rm(list=ls(all=TRUE))
x <- read.table("C:\\Users\\danyang\\Dropbox\\588\\HTF\\data\\train4.data",sep=",")
x <- as.matrix(x)
y <- array(x,c(nrow(x),16,16))
y <- y[,,16:1]
x1 <- matrix(y,ncol=256)
windows()
par(mfrow = c(4,4),mar=rep(1,4))
for(i in 1:16){
  image(y[i,,],col=gray((32:0)/32))
}
# load 9
x <- read.table("C:\\Users\\danyang\\Dropbox\\588\\HTF\\data\\train9.data",sep=",")
x <- as.matrix(x)
y <- array(x,c(nrow(x),16,16))
y <- y[,,16:1]
x2 <- matrix(y,ncol=256)
windows()
par(mfrow = c(4,4),mar=rep(1,4))
for(i in 1:16){
  image(y[i,,],col=gray((32:0)/32))
}
# combine
x <- rbind(x1,x2)
x.mean <- apply(x,2,mean)
# center
x.demean <- x-rep(1,nrow(x))%*%t(x.mean)
x.svd <- svd(x.demean)

# mean image and 8 eigen-images
windows()
par(mfrow = c(3,3),mar=rep(1,4))
image(matrix(x.mean,16),col=gray((32:0)/32))
for(i in 1:8){
  image(matrix(-x.svd$v[,i],16),col=gray((32:0)/32))
}
# scores
plot(x.svd$u[,1],x.svd$u[,2],col = c(rep(2,nrow(x1)),rep(3,nrow(x2))))
# singular values
plot(x.svd$d)
# approximation
r <- 10
x.reconstruct <- rep(1,nrow(x)) %*% t(x.mean) + 
                 x.svd$u[,1:r] %*% diag(x.svd$d[1:r]) %*% t(x.svd$v[,1:r])
# images of reconstruction
windows()
par(mfrow = c(4,4),mar=rep(1,4))
for(i in 1:8){
  image(matrix(x.reconstruct[i,],16),col=gray((32:0)/32))
}
for(i in 1:8){
  image(matrix(x.reconstruct[i+nrow(x1),],16),col=gray((32:0)/32))
}
###########################################################################

###########################################################################
# Three Gaussian example
rm(list=ls())
p <- 128
n <- 300
# mean for three classes
theta1 <- c(3,0,rep(0,p-2))
theta2 <- c(0,-3,rep(0,p-2))
theta3 <- c(-1,0,rep(0,p-2))
# the mean matrix
mu <- rbind(matrix(rep(theta1,n/3),nrow=n/3,byrow=TRUE),
            matrix(rep(theta2,n/3),nrow=n/3,byrow=TRUE),
            matrix(rep(theta3,n/3),nrow=n/3,byrow=TRUE))
# observed data, with cov matrix = diag(.5,2,...,2)
x <- mu + matrix(rnorm(n*p),n) %*% diag(c(.5,2,rep(2,p-2)))

# oracle projection; the first two coordinates
pdf(file="LDA_oracle.pdf",height=8,width=8)
plot(x[,1],x[,2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)))
points(mu[,1],mu[,2],pch=19,cex=2)
dev.off()

# PCA by SVD
x.svd <- svd(x)
# PCA projection
pdf(file="LDA_PCA.pdf",height=8,width=8)
plot(x.svd$u[,1]*x.svd$d[1],x.svd$u[,2]*x.svd$d[2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)))
dev.off()
# PCA loadings
pdf(file="LDA_PCA_loading.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(x.svd$v[,1])
plot(x.svd$v[,2])
plot(x.svd$v[,3])
plot(x.svd$v[,4])
dev.off()
# PCA spectral
pdf(file="LDA_PCA_spectral.pdf",height=8,width=8)
plot(x.svd$d)
dev.off()

# compute mean estimate: overall and group means
mu.hat <- apply(x,2,mean)
mu.1.hat <- apply(x[1:(n/3),],2,mean)
mu.2.hat <- apply(x[(n/3+1):(2*n/3),],2,mean)
mu.3.hat <- apply(x[(2*n/3+1):n,],2,mean)

# between class covariance
S.b <- ((n/3)*(mu.1.hat-mu.hat)%*%t(mu.1.hat-mu.hat)+
       (n/3)*(mu.2.hat-mu.hat)%*%t(mu.2.hat-mu.hat)+
       (n/3)*(mu.3.hat-mu.hat)%*%t(mu.3.hat-mu.hat))/(n-1)
# within class covariance
S.w <- (t(x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) %*% (x[1:(n/3),] - rep(1,n/3)%*% t(mu.1.hat)) +
       t(x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) %*% (x[(n/3+1):(2*n/3),] - rep(1,n/3)%*% t(mu.2.hat)) +
       t(x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)) %*% (x[(2*n/3+1):n,] - rep(1,n/3)%*% t(mu.3.hat)))/(n-3)
# total variance
S.t <- t(x - rep(1,n)%*% t(mu.hat)) %*% (x - rep(1,n)%*% t(mu.hat))
# relation 
# S.t - S.b*(n-1) - S.w*(n-K) = 0
max(abs(S.t - S.b*(n-1) - S.w*(n-3)))

# define relative matrix
S <- solve(S.w) %*% S.b
# eigen decomp. of S
S.eig <- svd(S)

# FDA projection
pdf(file="LDA_LDA.pdf",height=8,width=8)
plot(x %*% S.eig$v[,1], x%*% S.eig$v[,2], col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)))
dev.off()

# FDA loadings
pdf(file="LDA_LDA_loading.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(S.eig$v[,1])
plot(S.eig$v[,2])
plot(S.eig$v[,3])
plot(S.eig$v[,4])
dev.off()

# FDA spectral
pdf(file="LDA_LDA_spectral.pdf",height=8,width=8)
plot(S.eig$d)
dev.off()

# three projections together
pdf(file="LDA.pdf",height=3,width=8)
par(mfrow=c(1,3))
plot(x[,1],x[,2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)),main="oracle projection")
points(mu[,1],mu[,2],pch=19,cex=2)
plot(x.svd$u[,1]*x.svd$d[1],x.svd$u[,2]*x.svd$d[2],col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)),main="PCA projection")
plot(x %*% S.eig$v[,1], x%*% S.eig$v[,2], col=rep(2:4,rep(n/3,3)),pch=rep(2:4,rep(n/3,3)),main = "LDA projection")
dev.off()



