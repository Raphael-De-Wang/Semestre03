#!env Rscript

library(mclust)

loadData <- function(fname="obeses_temoins.txt") {
    X <- as.matrix(t(read.table(fname)))
    Y <- c(rep(1,25),rep(2,10))
    return(list("X"=X,"Y"=Y))
}

standardEM <- function(X,Y) {
    msEst <- mstep(modelName = "VII", data = X[,1:100], z = unmap(Y))
    em(modelName=msEst$modelName, data= X[,1:100], parameters=msEst$parameters)
    do.call("em", c(list(data = X), msEst))
}

data <- loadData()
standardEM(data$X,data$Y)


