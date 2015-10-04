#!env Rscript

pdf(file=ifelse(FALSE, "tp02_Q1Q2Q3_plots_summary.pdf", "tp02_Q1Q2Q3_plots_summary.pdf"))
attach(mtcars,warn.conflicts = FALSE)
# layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE))

#### #### Q1 #### ####

# import cluster package
library(cluster)

# set random seed
set.seed(19)

# artifical random data [x,y], size of 100
data <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))

# set col names
colnames(data) <- c("x", "y")

# appliquer kmeans 
cl1 <- kmeans(data, 2)

# visulization
plot(data, col = cl1$cluster)
points(cl1$centers, col = 1:2, pch = 8, cex=2)

#### #### Q2 #### ####

# Hierarchical cluster analysis 
hc <- hclust(dist(data))
plot(hc)

# cutree
# k: an integer scalar or vector with the desired number of groups
# h: numeric scalar or vector with heights where the tree should be cut.
plot(hc,cutree(hc, k = 2))

#### #### Q3 #### ####

silVec <- c()
for ( i in 2 : (nrow(data)-1)) {
    cl <- cutree(hc, k=i)
    sil <- silhouette(cl,dist(data))
    silMean <- mean(sil[,3])
    silVec <- c(silVec,silMean)
}

plot(silVec)

plot(silhouette(cutree(hc, k = 2), dist(data)))
plot(silhouette(cutree(hc, k = 48), dist(data)))

dev.off()

#### Q4. Clustering de profils d'expression transcriptionnelle ####

pdf(file=ifelse(FALSE, "tp02_Q4_plots_summary.pdf", "tp02_Q4_plots_summary.pdf"))
attach(mtcars,warn.conflicts = FALSE)

library(Hmisc)
library(FunNet)
library(sna)
# https://cran.r-project.org/web/packages/FactoMineR/FactoMineR.pdf
library(FactoMineR)

data(obese)
upm <- apply(as.matrix(up.frame), 1, as.numeric)
dnm <- apply(as.matrix(down.frame), 1, as.numeric)
ud <- cbind(upm,dnm)[-1,]

rownames(ud) <- colnames(up.frame[,-1])
colnames(ud) <- paste("g",as.character(1:ncol(ud)),sep="")

indnames <- rownames(ud)
lnames <- indnames[1:11]
tnames <- indnames[12:length(indnames)]

#### Q4.1 ####
# PCA(ud[lnames,])
# PCA(ud[tnames,])
res.pca <- PCA(ud,graph = F)
# selon le plot, on voit trois groupe de points de nuage. En haut a gauche, les points presentent les echantillons qui ne sont pas obeses. 
plot.PCA(res.pca)

#### Q4.2 ####

# Hierarchical cluster analysis 
res.hcpc.3Dmap <- HCPC(res.pca,graph=F)
plot(res.hcpc.3Dmap)

res.hcpc.map <- HCPC(res.pca,graph=F)
plot(res.hcpc.map,choice="map")

hc <- hclust(dist(ud))
plot(hc)

# cutree
plot(hc,cutree(hc, k = 3))

# silhouette
silVec <- c()
for ( i in 2 : (nrow(ud)-1)) {
    cl <- cutree(hc, k=i)
    sil <- silhouette(cl,dist(ud))
    silMean <- mean(sil[,3])
    silVec <- c(silVec,silMean)
}

plot(silVec,main="Silhouette Vectors")
ct <- cutree(hc, k = 3)
plot(silhouette(ct, dist(ud)))

plot.PCA(res.pca, choix="ind", select=names(ct[ct==1]), col.ind="blue",new.plot = F)
plot.PCA(res.pca, choix="ind", select=names(ct[ct==2]), col.ind="red",new.plot = F)
plot.PCA(res.pca, choix="ind", select=names(ct[ct==3]), col.ind="green",new.plot = F)

#### Q4.3 kmeans ####
cl2 <- kmeans(ud, 3)
# print(cl2$cluster)
plot(ud, col = cl2$cluster)
points(cl2$centers, col = 1:2, pch = 8, cex=2)

dev.off()
