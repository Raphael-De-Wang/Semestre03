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

glm.binomial <- Q1.SimulezJeu()

#### La regression logistique binaire

randInitTheta <- function (numFeature=2) {
    return(as.matrix(rnorm(numFeature+1),row=1))
}

vectorXmatrix <- function (theta,X) {
    return(unlist(lapply(c(1:nrow(X)),function(n,theta,X){t(theta)%*%X[n,]},theta,X)))
}

logVrai <- function (d, theta) {
    # remove y
    x <- as.matrix(d$frame)[,-1]
    y <- as.matrix(d$y)
    thetaX <- vectorXmatrix(theta,x)
    return(sum(log(exp(thetaX)+1)-y*thetaX))
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
    print(vectorXmatrix((y-p(theta,x)),t(x)))
}

computeHessianMtx <- function (theta,x,y) {
    mxv <- function(i,m,y){return(m[i,]*y)}
    print(matrix(unlist(lapply(1:nrow(x),mxv,tcrossprod(x),p(theta,x))),nrow=nrow(x))) ####
    # unlist(lapply(1:nrow(x),mxv,tcrossprod(x),p(theta,x)))
}

basefuns.testcase <- function(){
    # build fake data 
    fdata <- fakeData(4)
    plot(fdata$x1,fdata$x2)

    # random model
    theta <- randInitTheta()
    print(theta)

    # log vraisemblance
    lv <- logVrai(fdata,theta)
    print(sprintf("log vraisemblance : %f",lv))

    x <- matrix(unlist(fdata$frame[,-1]), nrow=fdata$numEchant)
    y <- matrix(unlist(fdata$y))
    # probability
    pRes <- p(theta,x)
    print(fdata$y)
    print(sprintf("y=1 probability : %f", pRes))

    # Compute the first derivative
    firstDerivative(theta,x,y)

    # Compute the Hessian matrix
    computeHessianMtx(theta,x,y)
}

NewtonRaphon <- function (d,eps=1.e-5,maxIter=500) {
    theta <- randInitTheta()
    lvlist <- c()
    # iterate to the convergence
    while ( (length(lvlist)<2 || abs(lvlist[1]-lvlist[2])>eps) && (length(lvlist) < maxIter) ) {
        print(sprintf("iteration number [%d] [%f]", length(lvlist), lvlist[1]))  
        lv <- logVrai(d,theta)
    }
}

basefuns.testcase()

# NewtonRaphon(fakeData(4))
