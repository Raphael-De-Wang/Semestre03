#!env python

library(numDeriv)

fakeData <- function (n=1000) {
    set.seed(666)
    x0 = rep(1,n)
    x1 = rnorm(n)
    x2 = rnorm(n)
    z = 1 + 2*x1 + 3*x2
    pr = 1/(1+exp(-z))
    y = rbinom(n,1,pr)
    # convetion : col 0 is y
    df = data.frame(y=y,x0=x0,x1=x1,x2=x2)
    return(list(x0=x0,x1=x1,x2=x2,y=y,z=z,pr=pr,frame=df,d=2,numEchant=n,nomFeauture=2))
}

fakeData2 <- function () {
    x0 = rep(1,4)
    x1 = c(1,2,-3,-4)
    x2 = c(1,-2,3,-4)
    y  = c(0,0,1,1)
    df = data.frame(y=y,x0=x0,x1=x1,x2=x2)
    return(list(x0=x0,x1=x1,x2=x2,y=y,frame=df,numEchant=4,nomFeauture=2))
}

fakeData3 <- function(n1,n2) {
    randOrder <- sample(1:(n1+n2))
    x0 = rep(1,n1+n2)
    x1 <- c(rnorm(n1, mean = -5),rnorm(n2, mean =  10))[randOrder]
    x2 <- c(rnorm(n1, mean = -5),rnorm(n2, mean =  10))[randOrder]
    y <- c(rep(0,n1),rep(1,n2))[randOrder]
    df = data.frame(y=y,x0=x0,x1=x1,x2=x2)
    numEchant <- n1 + n2
    numFeature <- 2
    return(list(x0=x0,x1=x1,x2=x2,y=y,frame=df,"numEchant"=numEchant,"numFeature"=numFeature))
}


Q1.SimulezJeu <- function () {
    pdf(file=ifelse(FALSE, "tp04_Q1SimulezJeu.pdf", "tp04_Q1SimulezJeu.pdf"))
    attach(mtcars,warn.conflicts = FALSE)
    fData <- fakeData(100)
    glm.binomial <- glm( fData$y~fData$x1+fData$x2,data=fData$frame,family="binomial")
    color <- fData$y
    color[color==1] <- 'blue'
    color[color==0] <- 'red'
    plot(x=fData$x1, y=fData$x2, pch = 3, type = "p", col=color)
    for (i in 1:6) {
        plot(glm.binomial,which=c(i))
    }
    dev.off()
    return(glm.binomial)
}

#### La regression logistique 

randInitTheta <- function (numFeature=2) {
    return(matrix(c(1,-1,1),ncol=1))
}

logVrai <- function (d, theta) {
    x <- as.matrix(d$frame)[,-1]
    y <- as.matrix(d$y)
    return(sum(unlist(lapply(1:nrow(x),function(n,theta,x,y){y[n]*t(theta)%*%x[n,]-log(1+exp(t(theta)%*%x[n,]))},theta,x,y))))
}

vectorXmatrix <- function (theta,X) {
    return(unlist(lapply(c(1:nrow(X)),function(n,theta,X){t(theta)%*%X[n,]},theta,X)))
}

p <- function (theta, x, y = 1) {
    expThetaX <- exp(vectorXmatrix(theta,x))
    if ( y == 1 ) {
        return(expThetaX/(1+expThetaX))
    } else {
        return(1/(1+expThetaX))
    }
}

firstDerivative <- function (theta,x,y) {
    a <- lapply(1:nrow(x),function(n,prob,x,y){return(x[n,]*(prob[n]-y[n]))},p(theta,x),x,y)
    return(colSums(matrix(unlist(a),nrow = nrow(x), byrow = TRUE)))
}

computeHessianMtx <- function (theta,x) {
    mxv <- function(n,x,prob) { return(tcrossprod(x[n,])*prob[n]*(1-prob[n])) }
    H <- matrix(rep(0,ncol(x)^2),nrow=ncol(x))
    for ( hssn in lapply(1:nrow(x),mxv,x,p(theta,x))) {
        H <- H + hssn
    }
    return(H)
}

dvp <- function (d,theta) {
    theta <- theta/theta[2]
    plot(d$x1,d$x2,type="n")
    text(d$x1,d$x2,labels = d$y,col = c('red','blue','green')[d$y+1])
    abline(a = -theta[1]/theta[3], b = -theta[2]/theta[3] , col = 2)
}

basefuns.testcase <- function(n){
    # build fake data 
    # fdata <- fakeData(n)
    fdata <- fakeData2()
    
    # random model
    theta <- randInitTheta()
    print(theta)

    # plot data and theta
    dvp(fdata,theta)
    
    # log vraisemblance
    lv <- logVrai(fdata,theta)
    print(sprintf("log vraisemblance : %f",lv))

    x <- matrix(unlist(fdata$frame[,-1]), nrow=fdata$numEchant)
    y <- matrix(unlist(fdata$y))

    # probability
    pRes <- p(theta,x)
    print(fdata$y)
    print(sprintf("y=1 probability : %f", pRes))

    print("Compute the first derivative : ")
    print(firstDerivative(theta,x,y))

    print("Compute the Hessian matrix : ")
    print(computeHessianMtx(theta,x))

}


matrixXvector <- function(m,v) {
    return(matrix(unlist(lapply(1:nrow(m),function(n,m,v){m[n,]%*%v},m,v))))
}


NewtonRaphon <- function (d,eps=1.e-5,maxIter=500) {
    lvlist <- c()
    theta <- randInitTheta()
    x <- matrix(unlist(d$frame[,-1]), nrow=d$numEchant)
    y <- matrix(unlist(d$y))
    # iterate to the convergence
    while ( (length(lvlist) < 2 || abs(lvlist[1]-lvlist[2])>eps) && (length(lvlist) < maxIter) ) {
        dvp(d,theta)
        lvlist <- c(logVrai(d,theta),lvlist)
        print(sprintf("iteration number [%d] [%f]", length(lvlist), lvlist[1]))
        print(theta)
        sh <- solve(computeHessianMtx(theta,x))
        fd <- firstDerivative(theta,x,y)
        theta <- theta - matrixXvector(sh,fd)
    }
    return(theta)
}

NewtonRaphon.testcase <- function(n) {
    pdf(file=ifelse(FALSE, "Newtown_test.pdf", "Newtown_test.pdf"))
    attach(mtcars,warn.conflicts = FALSE)
    d <- fakeData(n)
    # d <- fakeData2()
    # d <- fakeData3(10,20)
    theta <- NewtonRaphon(d)
    dvp(d,theta)
    dev.off()
}

# glm.binomial <- Q1.SimulezJeu()
# basefuns.testcase(40)
NewtonRaphon.testcase(50)
