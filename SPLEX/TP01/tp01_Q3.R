#!env Rscript
library(matrixStats)
library(impute)
library(samr) 
library(corrplot)

pdf("Q3_summary.pdf")
attach(mtcars,warn.conflicts = FALSE)
layout(matrix(c(1,2,5,3,3,4), 2, 3, byrow = TRUE))

# load data
x <- read.delim("obeses_temoins.txt")

# header list
figList <- colnames(x)
obeses  <- figList[1:25]
temoins <- figList[26:35]
geneid  <- row.names(x)
# correlation
corrplot(cor(x[obeses]),type="upper", order="hclust", tl.col="black", tl.srt=45, method="ellipse")
corrplot(cor(x[temoins]),type="upper", order="hclust", tl.col="black", tl.srt=45, method="pie")

# samr
delta = 0.6756
data = list(x=as.matrix(x),y=c(rep(1,25),rep(2,10)),geneid=geneid, genenames=geneid, logged2=TRUE)
samr.obj <- samr(data, resp.type="Two class unpaired", nperms=300)
samr.plot(samr.obj,delta)
delta.table <- samr.compute.delta.table(samr.obj)
siggenes.table <- samr.compute.siggenes.table(samr.obj, delta, data, delta.table)

write.table(siggenes.table$genes.up,'obeses_temoins_sig_up.txt',sep='\t',quote=F)
write.table(siggenes.table$genes.lo,'obeses_temoins_sig_lo.txt',sep='\t',quote=F)

pv <- c()
for ( i in 1:nrow(x) ) {
    pv <- c(pv, t.test(x[i,1:10],x[i,11:20])$p.value)
}

mlist <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
adjust <- new.env()
for ( m in mlist ) {
    adjust[[m]] <- p.adjust(pv, m)
    print(m)
    print(summary(adjust[[m]]))
}

boxplot(adjust[["holm"]],adjust[["hochberg"]],adjust[["hommel"]],adjust[["bonferroni"]],adjust[["BH"]],adjust[["BY"]],adjust[["fdr"]],adjust[["none"]],col = "lightgray", xaxt = "n",  xlab = "")
axis(1, labels = FALSE)
text(x = seq_along(mlist), labels = mlist, y = par("usr")[3], srt = 45, adj = 1, xpd = TRUE)

mtx  <- matrix(nrow=8,ncol=nrow(x))
for ( m in 1:length(mlist) ) {
    mtx[m,1:nrow(x)]  <- adjust[[mlist[m]]]
}

rownames(mtx) <- mlist
corrplot(cor(t(mtx)),type="upper", order="hclust", tl.col="black", tl.srt=45, method="pie")

sseuil <- 0.05
# print(sum((adjust[["BY"]] < sseuil) & (adjust[["bonferroni"]] < sseuil) * 1))
byIDS   <- adjust[["BY"]] < sseuil
bonfIDS <- adjust[["bonferroni"]] < sseuil
commIDS <- byIDS & bonfIDS
byOnlyIDS <- byIDS & !commIDS
bonfOnlyIDS <- bonfIDS & !commIDS

print(summary(delta.table))

print("common genes (BY and bonferroni): ")
print(data$genenames[commIDS])
print("BY method only genes : ")
print(data$genenames[byOnlyIDS])
print("bonferroni method only genes : ")
print(data$genenames[bonfOnlyIDS])

print("Correlation of BY and bonferroni method : ")
print(cor(adjust[["BY"]],adjust[["bonferroni"]]))

# print(colnames(siggenes.table$genes.up))
samrGIDs <- c(siggenes.table$genes.lo[,2],siggenes.table$genes.up[,2])
byGIDs   <- data$genenames[byIDS]
bonfGIDs <- data$genenames[bonfIDS]

print("common genes (samr and BY): ")
print(intersect(samrGIDs,byGIDs))
print("common genes (samr and bonferroni): ")
print(intersect(samrGIDs,bonfGIDs))

dev.off()

