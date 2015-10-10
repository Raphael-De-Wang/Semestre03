#!env Rscript

library(mvtnorm)

loadData <- function(fname="obeses_temoins.txt") {
    X <- as.matrix(t(read.table(fname)))
    Y <- c(rep(1,25),rep(2,10))
    return(list("X"=X,"Y"=Y))
}

dmvN <- function(x,mu,var,log=F){
    inner <- function(v1,v2){
        return(sum(v1*v2))
    }
    normalization <- function(x,var) {
        d <- ncol(x)
        return(((2*pi)^(d/2))*sqrt(det(var)))
    }
    exponent <- function(x,mu,var){
        return(inner(crossprod((-0.5)*t(x-mu),solve(var)),(x-mu)))
    }
    if (!is.matrix(var) || !is.matrix(mu) || !is.matrix(x) ) {
        stop("[constraint of variable] x: matrix, mu: matrix, var: matrix. Check dnorm for Linear Gaussien.")}
    if (ncol(var)!=nrow(var)) {
        stop("var is not symmetric.")
    }
    if (log == FALSE) {
        return(exp(exponent(x,mu,var))/normalization(x,var))
    } else {
        return(exponent(x,mu,var) - log(normalization(x,var)))
    }
}

dmvNTestCase <- function(n=3,numiter=50,eps=1.e-12,log=F) {
    for (i in 1:numiter) {
        x <- matrix(runif(n),nrow=1)
        mu <- matrix(runif(n),nrow=1)
        var <- matrix(runif(n^2),nrow=n) * diag(n)
        if (abs(dmvnorm(x,mu,var,log) - dmvN(x,mu,var,log))>eps) {
            print(abs(dmvnorm(x,mu,var,log) - dmvN(x,mu,var,log)))
            stop("dmvN meet a error.")
        }
    }
    print("Test Case End without error.")
}

data <- loadData()
# dmvNTestCase(n=100,log=T)

LogVraisemblance <- function(X,alpha,mu,var) {
    sapply(1:nrow(X),dmvnorm,X,mu,var)
}
