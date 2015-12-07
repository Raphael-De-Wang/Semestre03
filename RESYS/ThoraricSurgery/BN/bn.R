#!env R

library(foreign)

library(infotheo)
library(bnlearn)
library(pcalg)

ts <- read.arff("../data/ThoraricSurgery.arff")

#### feature engineering ####

# numeric
for ( name in attributes(ts)$names ) {
    ts[,name] <- as.numeric(ts[,name])
}

# X <- ts[,-which( colnames(ts)=="Risk1Yr" )]
# Y <- ts$Risk1Yr

# discretization


gs.res <- gs(ts)
# summary(gs.res)
plot(gs.res)
# plot(hc(ts))

