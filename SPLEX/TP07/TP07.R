#env Rscript

library(MASS)

pdf(file=ifelse(FALSE, "tp07.pdf", "tp07.pdf"))
attach(mtcars,warn.conflicts = FALSE)

trset1 <- read.csv("data1_train.txt",header=T)
ttset1 <- read.csv("data1_test.txt",header=T)

trY <- read.csv("class_train.txt",header=F)
ttY <- read.csv("class_test.txt",header=F)

plot.tp07 <- function(dataX,Y) {
    mycol = c("red","blue")
    family = as.factor(t(Y))
    palette(mycol)
    plot(x=dataX[,1],y=dataX[,2],col=family,pch=3)
}

plot.tp07(trset1,trY)
plot.tp07(ttset1,ttY)

tauxErr = c()

####
trainSet <- data.frame(trset1,trY)
colnames(trainSet) = c("v1","v2","y")
testSet <- data.frame(ttset1)
colnames(testSet) = c("v1","v2")
z <- lda(y ~ ., trainSet, prior = c(1,1)/2, subset = 1:nrow(trainSet))
y <- predict(z, testSet)$class
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))


####
z <- qda(trset1, as.factor(t(trY)))
y <- predict(z,ttset1)$class
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))


####
library(class)
y <- knn(trset1, ttset1, t(trY), k = 3, l = 0, prob = FALSE, use.all = TRUE)
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))

#### 
library(e1071)
trainSet <- data.frame(trset1,trY)
colnames(trainSet) = c("v1","v2","y")
mGauss <- svm(y~., data=trainSet, gamma = 0.01)
y <- predict (mGauss, ttset1)
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))

## 
mLinear <- svm(y~., data=trainSet, gamma = 0.01, kernel="linear")
y <- predict (mLinear, ttset1)
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))

## 
mPoly <- svm(y~., data=trainSet, gamma = 0.01, kernel="polynomial")
y <- predict (mPoly, ttset1)
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))


####
library(rpart)
fit <- rpart(y~v1+v2, data=trainSet)
plot(fit)
text(fit, use.n=TRUE, all=TRUE, cex=.8)
plotcp(fit)
y <- predict(fit, testSet, type="class")
plot.tp07(ttset1,t(y))
tauxErr = c(tauxErr,sum(y != as.factor(t(ttY))) / nrow(ttY))


####
library(adabag)
ada <- function(trainSet,testSet,trY,ttY) {
    tauxErr = c()
    for ( niter in 1:10 ) {    
        adaboost <- boosting(y~., data=trainSet, boos=TRUE, mfinal=niter)
        adaboost.pred <- predict.boosting(adaboost, testSet)
        tauxErr <-c(tauxErr,niter =  sum(adaboost.pred$class != as.factor(t(ttY))) / nrow(ttY))
    }
    return(tauxErr)
}

res = replicate(5,ada(trainSet,testSet,trY,ttY))

print(res)

dev.off()
