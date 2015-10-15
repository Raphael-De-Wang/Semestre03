#!env Rscript

pdf(file=ifelse(FALSE, "tp03_Analyse_summary.pdf", "tp03_Analyse_summary.pdf"))
attach(mtcars,warn.conflicts = FALSE)

library(mclust)

loadData <- function(fname="obeses_temoins.txt") {
    X <- as.matrix(t(read.table(fname)))
    Y <- c(rep(1,25),rep(2,10))
    return(list("X"=X,"Y"=Y))
}

EMMixModel <- function(X,Yp,Y,modelName) {
    msEst <- mstep(modelName = "VVV", data = X[,1:1000], z = unmap(Yp))
    resEM <- em(modelName=msEst$modelName, data= X[,1:1000], parameters=msEst$parameters)
    return(resEM)
}

visualizeData <- function(X) {
    library(FactoMineR)
    res.pca <- PCA(data$X,graph = F)
    plot.PCA(res.pca)
    res.hcpc.3Dmap <- HCPC(res.pca,graph=F)
    plot(res.hcpc.3Dmap)
}

ClusterAnalysis <- function (X) {
    library(gridExtra)
    res <- Mclust(X,modelNames=c("VII","VVV"))
    grid.table(res$BIC)
    plot(res$BIC,what="BIC")
    plot(res$BIC,what="classification")
    plot(res$BIC,what="uncertainty")
    plot(res$BIC,what="density")
    return(res)
}

ClusterAnalysisBIC <- function (X) {
    resBIC <- mclustBIC(X,modelNames=c("VII","VVV"))
    Xsum <- summary(resBIC,X)
    mclust2Dplot(X,classification=Xsum$classification,parameters=Xsum$parameters)
    return(resBIC)
}

ClusterAnalysisPrior <- function (X) {
    resBIC <- mclustBIC(X,modelNames=c("VII","VVV"),prior=priorControl())
    Xsum <- summary(resBIC,X)
    mclust2Dplot(X,classification=Xsum$classification,parameters=Xsum$parameters)
    return(resBIC)
}


data <- loadData()

print("visualizeData")
visualizeData(data$X)

print("ClusterAnalysis")
ClusterAnalysis(data$X)

print("ClusterAnalysisBIC")
ClusterAnalysisBIC(data$X)

# ClusterAnalysisPrior(data$X)

dev.off()
