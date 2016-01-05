#!env Rscript

pdf(file=ifelse(FALSE, "FunNet_Res.pdf", "FunNet_Res.pdf"))
attach(mtcars,warn.conflicts = FALSE)

# 1. Exemple d'annotation fonctionnelle avec le package FunNet:
#    http://bioinformatics.oxfordjournals.org/content/24/22/2636.full
#    http://www.funnet.ws/

library(FunNet)

load("obese.RData")

# 1.1 Créer deux variables contenant les noms des gènes (les Entrez GeneID) sur et sous exprimés en les transformant en data.frame. 

# 1.2 Ajouter à chacune de ces data.frames 3 colonnes supplémentaires contenant que des valeurs nulles. 

surGIDs <- data.frame(up.frame,S0=rep(0,nrow(up.frame)),S1=rep(0,nrow(up.frame)),S2=rep(0,nrow(up.frame)))
sousGIDs <- data.frame(down.frame,S0=rep(0,nrow(down.frame)),S1=rep(0,nrow(down.frame)),S2=rep(0,nrow(down.frame)))

# 1.3 Appeler la fonction FunNet() en changent les paramètres : 

FunNet(org="HS",up.frame=surGIDs,down.frame=sousGIDs,ref.list=NULL,
       two.lists=TRUE,restrict=FALSE,go.bp=TRUE,kegg=TRUE,
       discriminant=TRUE,annot.method="specificity",fdr=NA,
       build.annot.net=FALSE,estimate.th=FALSE,level=1,
       hard.th=NA,soft.th=NA,topological=FALSE,coexp.method="spearman",
       keep.rdata=TRUE,zip=FALSE)

dev.off()
