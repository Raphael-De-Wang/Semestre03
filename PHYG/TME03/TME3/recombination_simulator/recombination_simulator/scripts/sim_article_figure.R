# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

rm(list=ls())

while (length(dev.list()) > 0)
  dev.off()
while (sink.number() > 0)
  sink()

setwd("~/Documents/PhyChro/sim")

source("scripts/sim_functions.R")

for (i in 1:2) {
  rlcLoadOnly(c("results_yeast1/benchmark-rdata.rda", "results_vert1/benchmark-rdata.rda")[[i]],
              c("incorrectEdgeLengths","correctEdgeLengths","allDist"))
  pdf(c("summary/sim-yeast.pdf","summary/sim-vert.pdf")[[i]])
  
  X <- tabulate(allDist+1)
  #l <- paste(0:(length(X)-1), "errors")
  #l[[1]] <- "No errors"
  #l[[2]] <- "1 error"
  l <- 0:(length(X)-1)
  z <- barplot(X,
               ylim=c(0,100),
               names.arg=l,
               col=c("white",rep("black", length(X)-1)),
               #main=paste("Number of correct/incorrect splits in trees -", c("Yeast","Vertebrates")[[i]]),
               main=NA,
               cex.names=1.5,
               cex.axis=1.5)
  text(z, X, labels=round(100*X/sum(X)), pos=3, xpd=TRUE, cex=1.5)
  
  rect(usrFromRelativeX(0.40), usrFromRelativeY(0.33), usrFromRelativeX(1.02), usrFromRelativeY(1), xpd=TRUE)
  par(fig=c(0.46, 1, 0.4, 0.9), new=TRUE, cex=0.8, mai=c(0.7,0.7,0.6,0.4), xpd=TRUE)
  
  X <- c(incorrectEdgeLengths,correctEdgeLengths)
  #breaks <- seq(0, bs+max(X, na.rm=TRUE), by=bs)
  breaks <- seq(0, max(X), length.out=50)
  h <- hist(X, breaks=breaks, plot=FALSE)
  h4 <- hist(correctEdgeLengths, breaks=breaks, plot=FALSE)
  plot(h, col="black", xlab="Number of events per branch", main="Branch lengths", ylab="Number of branches",
       cex.lab=1.5, cex.main=1.5)
  lines(h4, col="white")
#  legend("topright", fill=c("black","white"), legend=c("incorrect","correct"))
  
  dev.off()
}
