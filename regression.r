###########################################################################
#                           Name: Long T. Le
#                           Netid: longtle
###########################################################################
#Implement soft-thresholding function
sign <-function(x, lambda){
    y <- 0
    if (x > lambda) {
        y <- (x - lambda)
    }
    if (x < -lambda) {
        y <- x + lambda
    }
    return (y)
}

###########################################################################
#implement betabar function, calculate betaj based on existing beta
betabar <- function(beta, x, y, lambda,j) {
    n <- nrow(x)
    p <- ncol(x)
    firstpara <- 0
    for (i in 1:n){
        sum.xik.betak <- 0
        for (k in 1:p) {
            if (k != j) {
                sum.xik.betak <- (sum.xik.betak + x[i,k]*beta[k])
            }
        }
        firstpara <- firstpara + x[i,j]*(y[i] - sum.xik.betak)
    }
    return (sign(firstpara,lambda))
    
}

###########################################################################

# Coordinate Descent for LASSO in R (detail implementation)
mylasso <- function(x,y,lambda){
    x <- as.matrix(x)
    y <- as.vector(y)
    p <- ncol(x)
    tol <- 1e-5
    niter <- 100
    iiter <- 1
    dist <- 1
    beta.old <- vector(mode="numeric", length=p)
    
    beta.new <- beta.old
    while(iiter <= niter & dist > tol) {
        for(j in 1:p){
            beta.new[j] <- betabar(beta.new, x, y, lambda, j)
        }
        dist <- norm((beta.new - beta.old), "2")
        iiter <- iiter + 1
        beta.old <- beta.new
    }
    return(beta.new)
}

###########################################################################

#Sample example to test
x <- matrix(rnorm(10),5,2) # rnorm(): generate norm random variables
y <- rnorm(5)
x <- scale(x)
x <- x/norm(x[,1], "2")
y <- y - mean(y)
beta.hat <- mylasso(x,y,0.1)
beta.hat
beta.hat <- mylasso(x,y,1)
beta.hat
















