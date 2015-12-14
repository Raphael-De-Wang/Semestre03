#!env Rscript

library(bnlearn)
library(infotheo)
library(FactoMineR)

pdf(file=ifelse(FALSE, "tp09.pdf", "tp09.pdf"))

load('MetagenomicsExp.R')

# rs <- PCA(gene_ab,graph=FALSE,ncp=20)

# for ( i in 1:length(gene_ab[1,]) ) { 
#     if ( gene_ab[1,i] == rs$ind$coord[1,1] ) {
#         print(i)
#     } 
# }

# str(gene_ab)
# str(y)

# d <- discretize(gene_ab,nbins=5)
d <- data.frame(apply(gene_ab, 2, discretize, nbins=3))

for ( name in attributes(d)$names ) {
    d[,name] <- as.factor(d[,name])
}

hc.bn <- hc(d[,1:20])
plot(hc.bn)
title("Score-based structure learning algorithms : hill-climbing")
# print(modelstring(hc.bn))
plot(boot.strength(d[,1:20],algorithm="hc"))
# warnings()
# print(arc.strength(hc.bn, d[1:20]))

tabu.bn <- tabu(d[,1:20])
plot(tabu.bn)
title("Score-based structure learning algorithms : Tabu search (TABU)")
# print(modelstring(tabu.bn))

gs.bn <- gs(d[,1:20])
plot(gs.bn)
title("Constraint-based structure learning algorithms : Grow-Shrink (GS)")

iamb.bn <- iamb(d[,1:20])
plot(iamb.bn)
title("Constraint-based structure learning algorithms : the Incremental Association (IAMB)")

fast.iamb.bn <- fast.iamb(d[,1:20])
plot(fast.iamb.bn)
title("Constraint-based structure learning algorithms : the Fast Incremental Association (Fast-IAMB)")

inter.iamb.bn <- inter.iamb(d[,1:20])
plot(inter.iamb.bn)
title("Constraint-based structure learning algorithms : the Interleaved Incremental Association (Inter-IAMB) constraint-based algorithms.")

rsmax2.bn <- rsmax2(d[,1:20])
plot(rsmax2.bn)
title("Hybrid structure learning algorithms : Max-Min Hill Climbing (MMHC)")

mmhc.bn <- mmhc(d[,1:20])
plot(mmhc.bn)
title("More General 2-Phase Restricted Maximization (RSMAX2) Hybrid Algorithms")

# arc operations.
# set.arc(x, from, to, check.cycles = TRUE, debug = FALSE)
# drop.arc(x, from, to, debug = FALSE)
# reverse.arc(x, from, to, check.cycles = TRUE, debug = FALSE)

# edge (i.e. undirected arc) operations
# set.edge(x, from, to, check.cycles = TRUE, debug = FALSE)
# drop.edge(x, from, to, debug = FALSE)

# bn.boot

strength.res.boot <- boot.strength(d[,1:20],algorithm="hc")

plot(fast.iamb.bn)

for (ind in which(strength.res.boot$strength < 0.9)) {
    drop.arc(fast.iamb.bn, strength.res.boot$from[ind], strength.res.boot$to[ind])
    # print(c(strength.res.boot$from[ind], strength.res.boot$to[ind]))
}

plot(fast.iamb.bn)

dev.off()



