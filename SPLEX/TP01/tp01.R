#!env R

#### Q1 ####

# import samr library
library("samr") 

# init random number generator (RNG)
set.seed(100)   

# init artifical dataset, niveau d'expression, 20 patients, 1000 genes for each
x<-matrix(rnorm(1000*20),ncol=20)

# random index
dd<-sample(1:1000,size=100)

# noise
u<-matrix(2*rnorm(100),ncol=10,nrow=100)

# drip random noise
x[dd,11:20]<-x[dd,11:20]+u

# class label
y<-c(rep(1,10),rep(2,10))

# encapsulate data
data=list(x=x,y=y,geneid=as.character(1:nrow(x)),genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)

# Significance analysis of microarrays
# resp.type=c("Quantitative","Two class unpaired","Survival","Multiclass", "One class", "Two class paired",
#     "Two class unpaired timecourse", "One class timecourse", "Two class paired timecourse", "Pattern discovery"),
samr.obj<-samr(data, resp.type="Two class unpaired", nperms=100)

delta=.4

# plot
# samr.plot(samr.obj,delta)

# Computes tables of thresholds, cutpoints and corresponding False Discovery rates for SAM 
# delta.table <- samr.compute.delta.table(samr.obj)

# Computes significant genes table, starting with samr object "samr.obj" and delta.table "delta.table"
# siggenes.table<-samr.compute.siggenes.table(samr.obj, delta, data, delta.table)

#### Q2 ####
pv <- c()
for ( i in 1:nrow(x) ) {
    pv <- c(pv, t.test(x[i,1:10],x[i,11:20])$p.value)
}

mlist <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
adjustSum <- new.env()
for ( m in mlist ) {
    adjustSum[[m]] <- summary(p.adjust(pv, m))
}
