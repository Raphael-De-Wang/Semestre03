#!env Rscript

library(mvtnorm)

loadData <- function(fname="obeses_temoins.txt") {
    X <- as.matrix(t(read.table(fname)))
    X <- X[,1:2]
    Y <- c(rep(1,25),rep(2,10))
    numEchant <- nrow(X)
    numFeature <- ncol(X)
    return(list("X"=X,"Y"=Y,"numEchant"=numEchant,"numFeature"=numFeature))
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

dmvN.testcase <- function(n=3,numiter=50,eps=1.e-5,log=F) {
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

LogVraisemblance <- function(data,models) {
    alpha <- models$alpha
    mu <- models$mu
    var <- models$var
    z <- matrix(unlist(models$z),nrow = models$classNum, byrow = T)
    classNum <- models$classNum
    lv <- 0
    p <- matrix(rep(0,data$numEchant*2), ncol=data$numEchant)
    for ( indEchant in 1:data$numEchant) {
        for ( indClass in 1:classNum ) {
            x <- t(as.matrix(data$X[indEchant,]))
            p[indClass,indEchant] <- dmvnorm(x,mu[indClass,],matrix(unlist(var[indClass]),nrow=data$numFeature,byrow = TRUE),log=T)
            lv <- lv + z[indClass,indEchant]*(p[indClass,indEchant] + log(alpha[indClass]))
        }
    }
    return(list(lv=lv,p=p))
}

modelRandInit <- function(data,classNum=2) {
    alpha <- rep(1/classNum,classNum)
    mu <- array(runif(classNum*data$numFeature,min=-1,max=1),c(classNum,data$numFeature))
    var <- lapply(rep(1,classNum),diag,nrow=data$numFeature)
    z <- lapply(rep(1/classNum,times=classNum),rep,times=data$numEchant)
    return(list("alpha" = alpha, "mu" = mu, "var" = var, "classNum" = classNum, "z" = z))
}

estimation <- function(P) {
    p <- P/colSums(P)
    return(list(1 - ( p[2,] - p[1,] > 0 ) * 1,( p[2,] - p[1,] > 0 ) * 1))
}

maximisation <- function(data,models) {
    z <- matrix(unlist(models$z),nrow = models$classNum, byrow = T)
    models$alpha <- rowSums(z)/data$numEchant
    models$mu    <- crossprod(t(z),data$X)/rowSums(z)
    upVarm <- function(data,mu,z,m) {
        var <- matrix(rep(0,data$numFeature^2),nrow=data$numFeature)
        for ( i in 1:data$numEchant ) {
            if ( z[m,i] == 1 ) {
                var <- var + ((data$X[i,] - mu[m,]) %*% t(data$X[i,] - mu[m,]))
            }
        }
        return(var/sum(z[m,]))
    }
    models$var   <- lapply(1:models$classNum,upVarm,data=data,mu=models$mu,z=z)
    return(models)
}

estimationMaximisation <- function(data,eps=1e-2,maxIter=200){
    lvlist <- c()
    models <- modelRandInit(data)
    modelsfinal <- models
    while ( ( length(lvlist) < 2 || abs(lvlist[1]-lvlist[2]) > eps ) && (length(lvlist) < maxIter) ) {
        lv <- LogVraisemblance(data, models)
        models$z <- estimation(lv$p)
        models <- maximisation(data,models)
        lvlist <- c(lv$lv,lvlist)
        if ( length(lvlist) > 2 && lvlist[1] > max(lvlist[2:length(lvlist)]) ) {
            modelsfinal <- models            
        }
    }
    print(lvlist)
    print(which.max(lvlist))
    return(modelsfinal)
}

estimationMaximisation.testcase <- function(n1=25,n2=10){
    artificalData <- function(n1,n2) {
        d1 <- c(rnorm(n1, mean = -10),rnorm(n2, mean =  10))
        d2 <- c(rnorm(n1, mean = -10),rnorm(n2, mean =  10))
        randOrder <- sample(1:(n1+n2))
        X <- t(matrix(unlist(c(d1[randOrder],d2[randOrder])), nrow = 2, byrow = TRUE))
        numEchant <- nrow(X)
        numFeature <- ncol(X)
        return(list("X"=X,"numEchant"=numEchant,"numFeature"=numFeature))
    }
    data <- artificalData(n1,n2)
    models <- estimationMaximisation(data)
    plot(data$X)
    points(models$mu,col='red',pch = 21:25)
}

# data <- loadData()
# dmvN.testcase(n=100,log=T)
# models <- estimationMaximisation(data)
estimationMaximisation.testcase(n1=250,n2=100)
