#!env Rscript

library(mclust)

source("tp03_impl.R")

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

DensityEstimationAnalyse <- function (X) {
    resDens <- densityMclust(X)
    plot(resDens, what = "diagnostic")
    plot(resDens)
    plot(resDens, faithful, col = "grey", nlevels = 10, drawlabels = FALSE)
    plot(resDens, type = "image", col = topo.colors(50))
    plot(resDens, type = "persp", col = grey(0.8))
    return(resDens)
}

DiscriminantAnalyse <- function(X,Y) {
    x <- X[,1:1000]
    msEst <- mstep(modelName = "VVV", data = x, z = unmap(Y)) 
    emRes <- em(modelName = msEst$modelName, data = x, parameters = msEst$parameters)
    return(emRes)
}

mclustAnalyse <- function(){

    data <- loadData()

    # pdf(file=ifelse(FALSE, "tp03_Analyse_summary.pdf", "tp03_Analyse_summary.pdf"))
    # attach(mtcars,warn.conflicts = FALSE)

    print("visualizeData")
    visualizeData(data$X)
    
    print("ClusterAnalysis")
    ClusterAnalysis(data$X)

    print("ClusterAnalysisBIC")
    ClusterAnalysisBIC(data$X)
    
    # print("ClusterAnalysisPrior")
    # ClusterAnalysisPrior(data$X)

    print("DensityEstimationAnalyse")
    DensityEstimationAnalyse(data$X)
    dev.off()
}

# mclustAnalyse()
