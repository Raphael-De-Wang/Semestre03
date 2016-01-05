#!env Rscript

library(foreign)
library(infotheo)
library(bnlearn)
library(pcalg)
library(igraph)
library(minet)
library(discretization)
library(optparse)
library(corrplot)
library(Hmisc)

#  option_list = list(
#  make_option(c("-in"), action="store", default=NULL, type='character', dest = data,
#               help="data file name"),
#  make_option(c("--numeric"), action="store_true", default=FALSE, type='character', dest=numeric,
#              help="numerization of data and persist dictionary of factors"),
#  make_option(c("--discret"), action="store_true", default=FALSE, type='character', dest=discret,
#               help="discretization of data"),
#  make_option(c("--csv"), action="store", default=NULL, type='character', dest=csv,
#              help="persist data in csv"),
# )

# opt = parse_args(OptionParser(option_list=option_list))
# ts <- read.arff(opt$data)

arff_data <- read.arff("data/ThoraricSurgery.arff")

#### feature engineering ####

attributes_dict <- function () {
    DGN <- vector(mode="list", length=8)
    names(DGN) <- c("DGN1", "DGN2", "DGN3", "DGN4", "DGN5", "DGN6", "DGN8")
    DGN[[1]] <- 1;DGN[[2]] <- 2;DGN[[3]] <- 3;DGN[[4]] <- 4
    DGN[[5]] <- 5;DGN[[6]] <- 6;DGN[[7]] <- 7

    PRE6 <- vector(mode="list",length=3)
    names(PRE6) <- c("PRZ0", "PRZ1", "PRZ2")
    PRE6[[1]] <- 1; PRE6[[2]] <- 2; PRE6[[3]] <- 3;
    
    PRE14 <- vector(mode="list",length=4)
    names(PRE14) <- c("OC11", "OC12", "OC13", "OC14")
    PRE14[[1]] <- 1; PRE14[[2]] <- 2; PRE14[[3]] <- 3; PRE14[[4]] <- 4;

    LOGIC <- vector(mode="list",length=2)
    names(LOGIC) <- c("F","T")
    LOGIC[[1]] <- 1; LOGIC[[2]] <- 2

    PRE7 <- LOGIC
    PRE8 <- LOGIC
    PRE9 <- LOGIC
    PRE10 <- LOGIC
    PRE11 <- LOGIC
    PRE17 <- LOGIC
    PRE19 <- LOGIC
    PRE25 <- LOGIC
    PRE30 <- LOGIC
    PRE32 <- LOGIC
    Risk1Yr <- LOGIC

    return(list(DGN=DGN, PRE6=PRE6, PRE7=PRE7, PRE8=PRE8, PRE9 = PRE9, PRE10 = PRE10, PRE11 = PRE11, PRE14 = PRE14, PRE17 = PRE17, PRE19 = PRE19, PRE25 = PRE25, PRE30 = PRE30, PRE32 = PRE32, Risk1Yr = Risk1Yr))

}


no_discrete <- function (arff_data, name) {
    return(arff_data[,name])
}


numeric <- function (arff_data, discretize) {
    m <- matrix(rep(0,dim(arff_data)[[1]] * dim(arff_data)[[2]]), nrow=dim(arff_data))
    colnames(m) <- colnames(arff_data)
    data.attrs <- attributes_dict()
    for ( name in names(arff_data) ) {
        if ( name %in% names(data.attrs) ) {
            for ( elem_attr in names(data.attrs[[name]]) ) {
                m[which(arff_data[,name]==elem_attr),name] <- data.attrs[[name]][[elem_attr]]
            }
        } else {
            m[,name] <- discretize(arff_data,name)
        }
    }
    return(m)
}


persist <- function (m,fname) {
    write.table(m,file=fname,sep="\t",row.names=FALSE)
}


visualize <- function (m,filename) {
    pdf(file=ifelse(FALSE,filename,filename))
    d <- density(m[,"PRE4"])
    hist(m[,"PRE4"],prob=TRUE, col="grey",main="PRE4 - Forced vital capacity",ylim=c(0,max(d$y)))
    lines(d,col="blue", lwd=2)
    d <- density(m[,"PRE5"])
    hist(m[,"PRE5"],prob=TRUE, col="grey",main="PRE5 - Volume that has been exhaled at the end\n of the first second of forced expiration",ylim=c(0,max(d$y)))
    lines(d,col="blue", lwd=2)
    d <- density(m[,"AGE"])
    hist(m[,"AGE"],prob=TRUE, col="grey",main="AGE",ylim=c(0,max(d$y)))
    lines(d,col="blue", lwd=2)
    corrplot(rcorr(m)$r,method="ellipse")
    dev.off()
}


m <- numeric(arff_data,no_discrete)
visualize(m,"row_data.pdf")

cont_m <- matrix(rep(0,4*470),ncol=4)
colnames(cont_m) <- c("PRE4","PRE5","AGE","Risk1Yr")
cont_m[,"PRE4"] <- m[,"PRE4"]
cont_m[,"PRE5"] <- m[,"PRE5"]
cont_m[,"AGE"] <- m[,"AGE"]
cont_m[,"Risk1Yr"] <- m[,"Risk1Yr"]

ki2_res <- chi2(cont_m)

discrete_m <- m
discrete_m[,"PRE4"] <- ki2_res$Disc.data[,"PRE4"]
discrete_m[,"PRE5"] <- ki2_res$Disc.data[,"PRE5"]
discrete_m[,"AGE"]  <- ki2_res$Disc.data[,"AGE"]

visualize(discrete_m, "discrete_data.pdf")

persist(discrete_m,"/tmp/m.csv")
