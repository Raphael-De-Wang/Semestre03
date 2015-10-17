#!env Rscript

library(mvtnorm)

set.seed(1000)

loadData <- function(fname="obeses_temoins.txt") {
    X <- as.matrix(t(read.table(fname)))
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
        stop("[constraint of variable] x: matrix, mu: matrix, var: matrix. Check dnorm for Linear Gaussien.")
    }
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
            stop(sprintf("dmvN meet a error [%f].", abs(dmvnorm(x,mu,var,log) - dmvN(x,mu,var,log))))
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
    alpha <- rep(1,classNum)
    mu <- lapply(1:data$numFeature, function(n){return(runif(classNum)*max(data$X[,n]))})
    mu <- as.array(matrix(unlist(mu),nrow=classNum,byrow=TRUE))
    var <- lapply(rep(1,classNum),diag,nrow=data$numFeature)
    z <- lapply(1:classNum,function(n,partition){return((partition==n)*1)},partition=replicate(data$numEchant,sample(1:classNum,1)))
    return(list("alpha" = alpha, "mu" = mu, "var" = var, "classNum" = classNum, "z" = z))
}

estimation <- function(P) {
    p <- prop.table(P,2)
    return(list(1 - ( p[2,] - p[1,] > 0 ) * 1,( p[2,] - p[1,] > 0 ) * 1))
}

maximisation <- function(data,models) {
    z <- matrix(unlist(models$z),nrow = models$classNum, byrow = T)
    models$alpha <- rowSums(z)/data$numEchant
    models$mu    <- (crossprod(t(z),data$X))/rowSums(z)
    upVarm <- function(data,mu,z,oldvar,indClass) {
        var <- matrix(unlist(oldvar[indClass]),nrow=data$numFeature,byrow = TRUE) 
        for ( i in 1:data$numEchant ) {
            if ( z[indClass,i] == 1 ) {
                var <- var + ((data$X[i,] - mu[indClass,]) %*% t(data$X[i,] - mu[indClass,]))
            }
        }
        var <- var/sum(z[indClass,])
        return(var)
    }
    models$var <- lapply(1:models$classNum,upVarm,data=data,mu=models$mu,z=z,oldvar=models$var)
    return(models)
}

estimationMaximisation <- function(data,eps=1e-7,maxIter=500) {
    options(digits=20)
    lvlist <- c()
    models <- modelRandInit(data)
    while ( ( length(lvlist) < 2 || abs(lvlist[1]-lvlist[2]) > eps ) && (length(lvlist) < maxIter) ) {
        print(sprintf("iteration number [%d] [%f]", length(lvlist), lvlist[1]))
        lv <- LogVraisemblance(data, models)
        if (is.na(lv$lv)) {
            warning(sprintf("Data cannot be estmated into %d classes",models$classNum))
            break
        }
        models$z <- estimation(lv$p)
        models <- maximisation(data,models)
        lvlist <- c(lv$lv,lvlist)
        # print(models$mu)
        # print(models$var)
    }
    return(models)
}

estimationMaximisation.testcase <- function(n1=25,n2=10) {
    artificalData <- function(n1,n2) {
        d1 <- c(rnorm(n1, mean = -5),rnorm(n2, mean =  10))
        d2 <- c(rnorm(n1, mean = -5),rnorm(n2, mean =  10))
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

naiveTestSuite <- function () {
    data <- loadData()
    # test 1 : validation of multiple dimesional gaussien
    dmvN.testcase(n=100,log=T)

    # test 2 : EM test on biclass unbalanced data
    estimationMaximisation.testcase(n1=2,n2=20)
}

# naiveTestSuite()

# data <- loadData()
# models <- estimationMaximisation(data)
