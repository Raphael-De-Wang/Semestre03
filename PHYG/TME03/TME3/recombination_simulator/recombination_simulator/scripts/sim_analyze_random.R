# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

# Clear variables except for params
rm(list=setdiff(ls(), "params"))

source("scripts/sim_functions.R")

while (length(dev.list()) > 0)
  dev.off()
while (sink.number() > 0)
  sink()

if (!exists("params")) {
  stop("This script should be called after params have been defined.\n")
}
for (paramName in c("outDirectory")) {
  if (is.null(params[[paramName]])) {
    stop(paste("Missing parameter:", paramName))
  }
}
rm(paramName)

allSimOutDir <- paste(params$outDirectory, "/simulations", sep="")

sink(paste(params$outDirectory, "/benchmark-random.txt", sep=""), split=TRUE)

simIndex <- 1
t0 <- proc.time()[["elapsed"]]
incorrect <- 0
correct <- 0
allDist <- c()

correctEdgeLengths <- c()
incorrectEdgeLengths <- c()

correctConfidence <- c()
incorrectConfidence <- c()

allEdgeEvents <- numeric()
allEdgeEventsParams <- numeric()
# allEdgeEventsCorrect format:
# -l: branch is incorrect and of length l
# +l: branch is correct and of length l
# 0: impossible case
allEdgeEventsCorrect <- numeric()
smallBranchThreshold <- 5

while (TRUE) {
  simOutDir <- paste(allSimOutDir, "/", simIndex, sep="")
  if (!file.exists(simOutDir)) {
    break
  }

  correctTree <- read.tree(paste(simOutDir, "/tree-newick.txt", sep=""))
  correctTree$node.label <- 1:correctTree$Nnode
  correctTreeUnrooted <- unroot(correctTree)
  
  randomTree <- correctTreeUnrooted
  randomTree$tip.label <- sample(randomTree$tip.label, size=length(randomTree$tip.label), replace=FALSE)
  
  d <- dist.topo(correctTreeUnrooted, randomTree, method="PH85")
  bi <- prop.part(correctTreeUnrooted)
  z <- prop.clades(correctTreeUnrooted, randomTree)
  nErrors <- sum(z == 0)
  stopifnot(d == 2*nErrors)
  nSpecies <- Ntip(correctTreeUnrooted)
  edgeIndices <- match(nSpecies + (1:length(z)), correctTreeUnrooted$edge[,2])
  splitEdgeLengths <- correctTreeUnrooted$edge.length[edgeIndices]
  stopifnot(sum(is.na(splitEdgeLengths)) == 1)
  stopifnot(is.na(splitEdgeLengths[[1]]))
  
  incorrectEdgeLengths <- c(incorrectEdgeLengths, splitEdgeLengths[z == 0])
  correctEdgeLengths <- c(correctEdgeLengths, splitEdgeLengths[z == 1 & !is.na(splitEdgeLengths)])
  
  v <- prop.clades(randomTree, correctTreeUnrooted)
  v <- v[a2b(2, length(v))]
  confidence <- as.numeric(randomTree$node.label[a2b(2, length(randomTree$node.label))])
  incorrectConfidence <- c(incorrectConfidence, confidence[v == 0])
  correctConfidence <- c(correctConfidence, confidence[v == 1])
  
  #cat("*** Simulation ", simIndex, " -> incorrect splits = ", nErrors, "  ", (if (nErrors>0) "INCORRECT" else "OK"), "\n", sep="")
  cat("*** Simulation", simIndex)
  if (nErrors > 0) {
    cat("\n")
  } else {
    cat(" OK\n")
  }
  
  simDataFile <- paste(simOutDir, "/sim-R-data.RData", sep="")
  if (file.exists(simDataFile)) {
    simData <- new.env()
    load(simDataFile, envir=simData)
    M <- dist.nodes(simData$apeTreeE)
    for (i in a2b(2, length(z))) {
      ei <- edgeIndices[[i]]
      # Edge is a->b in unrooted tree
      a <- correctTreeUnrooted$edge[[ei,1]] - nSpecies
      b <- correctTreeUnrooted$edge[[ei,2]] - nSpecies
      l <- splitEdgeLengths[[i]]
      # Nodes ar and br are the respective node indices in rooted tree
      ar <- correctTreeUnrooted$node.label[[a]] + nSpecies
      br <- correctTreeUnrooted$node.label[[b]] + nSpecies
      # Distance between ar and br should be the same in both trees
      # otherwise there is a bug.
      stopifnot(M[[ar,br]] == l)
      rootedEdges <- which(simData$apeTreeE$edge[,1] == ar & simData$apeTreeE$edge[,2] == br)
      if (length(rootedEdges) == 0) {
        # Edge ei is in fact 2 edges in rooted tree
        rootedEdges <- which(simData$apeTreeE$edge[,2] == ar | simData$apeTreeE$edge[,2] == br)
        stopifnot(length(rootedEdges) == 2)
      }
      # rootedEdges now contains 1 or 2 edges in the rooted tree that correspond to the split
      stopifnot(sum(simData$apeTreeE$edge.length[rootedEdges]) == l)
      
      edgeEvents <- unlist(simData$events[rootedEdges])
      edegeEventsParams <- unlist(simData$eventsParams[rootedEdges])
      
      allEdgeEvents <- c(allEdgeEvents, edgeEvents)
      allEdgeEventsParams <- c(allEdgeEventsParams, edegeEventsParams)
      
      if (z[[i]] == 0) {
        cat("  Incorrect branch with length", l, "for split ")
        cat("[", paste(attr(bi, "labels")[bi[[i]]], collapse=","), "] vs others\n", sep="")
        cat("    with events:", paste(eventsAsChar(edgeEvents, edegeEventsParams), collapse=" "), "\n")
        allEdgeEventsCorrect <- c(allEdgeEventsCorrect, rep(-l, length(edgeEvents)))
      } else if (z[[i]] == 1) {
        allEdgeEventsCorrect <- c(allEdgeEventsCorrect, rep(l, length(edgeEvents)))
      }
    }
  }

  
  if (d > 0) {
    incorrect <- incorrect + 1
  } else {
    correct <- correct + 1
  }
  allDist <- c(allDist, nErrors)
  
  simIndex <- simIndex + 1
}

if (simIndex == 1) {
  cat("nothing to analyze\n")
} else {
  
  N <- correct + incorrect
  
  cat("Correct trees: ", correct, "/", N, " = ", 100*(correct/N), "%\n")
  
  cat("Details:\n")
  r <- range(allDist)
  for (i in r[[1]]:r[[2]]) {
    k <- sum(allDist == i)
    cat(" - ", k, " tree", (if (k>1) "s" else ""), " with ", i, " incorrect split", (if (i>1) "s" else ""), "\n", sep="")
  }
  
  cat("Total:", N, "trees\n")
  
  pdf(paste(params$outDirectory, "/plot-benchmark-random.pdf", sep=""), width=11, height=8)
  
  X <- tabulate(allDist+1)
  barplot(X,
          names.arg=paste(0:(length(X)-1), "\n", paste(round(100*X/sum(X)), "%"), sep=""),
          col=c("forestgreen",rep("orange", length(X)-1)),
          main=paste(paste(params$outDirectory, " - Number of incorrect splits in", length(allDist), "trees")),
          cex.names=0.9)
  
  X <- c(incorrectEdgeLengths,correctEdgeLengths)
  C <- c(rep(1, length(incorrectEdgeLengths)), rep(2, length(correctEdgeLengths)))
  layout(matrix(1:2, ncol=2))
  
  breaks <- seq(0, 5+max(X, na.rm=TRUE), by=5)
  h <- hist(X, breaks=breaks, plot=FALSE)
  
  h2 <- hist(incorrectEdgeLengths, breaks=breaks, plot=FALSE)
  
  if (length(incorrectEdgeLengths) > 0) {
    plot(h2, col="orange", main="Incorrect edges lengths", xlim=range(incorrectEdgeLengths), xlab="Number of events on branch", ylab="Number of branches")
  }
  
  h3 <- hist(correctEdgeLengths, breaks=breaks, plot=FALSE)
  
  #plot(h3, col="cyan", main="correct edges lengths")
  
  #
  # rlcCompareHist(X,
  #               c(rep(1,length(incorrectEdgeLengths)),rep(2,length(correctEdgeLengths))),
  #               breaks=breaks,
  #               col=c("orange","cyan"),
  #               classes=c("incorrect splits","correct splits"),
  #               relative=TRUE)
  
  h2$counts <- h2$counts / h$counts
  h2$counts[h$counts == 0] <- 0
  h3$counts <- h3$counts / h$counts
  h3$counts[h$counts == 0] <- 0
  
  h1 <- h2
  h1$counts <- h2$counts + h3$counts
  if (length(incorrectEdgeLengths) > 0) {
    plot(h1, col="orange", xlim=range(incorrectEdgeLengths), xlab="Number of events on branch", main=params$outDirectory)
    lines(h3, col="forestgreen")
  }
  
  layout(1)
  
  h4 <- hist(correctEdgeLengths, breaks=breaks, plot=FALSE)
  plot(h, col="orange", xlab="Number of events on branch", main=paste(params$outDirectory, "- Branch lengths"), ylab="Number of branches")
  lines(h4, col="forestgreen")
  legend("topright", fill=c("forestgreen", "orange"), legend=c("correct split", "incorrect split"))
  
  allConfidence <- c(correctConfidence, incorrectConfidence)
  hc <- hist(allConfidence, breaks=20, plot=FALSE)
  hci <- hist(incorrectConfidence, breaks=hc$breaks, plot=FALSE)
  hcc <- hist(correctConfidence, breaks=hc$breaks, plot=FALSE)
  plot(hc, col="orange", xlab="Confidence score", main=paste(params$outDirectory, "- PhyChro tree - branch confidence scores"), ylab="Number of branches")
  lines(hcc, col="forestgreen")
  legend("topleft", fill=c("forestgreen", "orange"), legend=c("correct split", "incorrect split"))
  
  plot(hci, col="orange", xlab="Confidence score", main=paste(params$outDirectory, "- PhyChro tree - branch confidence scores"), ylab="Number of branches")
  legend("topright", fill=c("orange"), legend=c("incorrect split"))
  
  if (length(allEdgeEvents) > 0) {
    
    mainTitle <- "Event types"
    correctEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect > 0]
    incorrectEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect < 0]
    M <- matrix(c(tabulate(incorrectEdgeEvents, nbins=7), tabulate(correctEdgeEvents, nbins=7)), nrow=2, ncol=7, byrow=TRUE)
    MC <- M
    M[1,] <- M[1,]/length(incorrectEdgeEvents)
    M[2,] <- M[2,]/length(correctEdgeEvents)
    z <- barplot(M, beside=TRUE, col=c("orange", "forestgreen"), names.arg=eventNames)
    legend("topright", fill=c("forestgreen", "orange"),
           legend=paste(c("correct splits", "incorrect splits"), " (", c(length(correctEdgeEvents), length(incorrectEdgeEvents)), " events)", sep=""))
    title(main=mainTitle)
    PV <- c()
    for (j in 1:7) {
      MF <- cbind(MC[,j], rowSums(MC[,-j]))
      fisherResult <- fisher.test(MF, alternative="two.sided")
      cat("Fisher test", eventNames[[j]], "vs branch correctness p-value =", fisherResult$p.value, "\n")
      PV <- c(PV, fisherResult$p.value)
    }
    text(x=colMeans(z), y=usrFromRelativeY(-0.1), xpd=TRUE, labels=paste0("p=", round(PV, 2)))

    
    correctEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect > 0 & abs(allEdgeEventsCorrect) <= smallBranchThreshold]
    incorrectEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect < 0 & abs(allEdgeEventsCorrect) <= smallBranchThreshold]
    M <- matrix(c(tabulate(incorrectEdgeEvents, nbins=7), tabulate(correctEdgeEvents, nbins=7)), nrow=2, ncol=7, byrow=TRUE)
    MC <- M
    if (sum(M) > 0) {
      M[1,] <- M[1,]/length(incorrectEdgeEvents)
      M[2,] <- M[2,]/length(correctEdgeEvents)
      z <- barplot(M, beside=TRUE, col=c("orange", "forestgreen"), names.arg=eventNames)
      legend("topright", fill=c("forestgreen", "orange"),
             legend=paste(c("correct splits with <= ", "incorrect splits with <= "), smallBranchThreshold, " events", " (", c(length(correctEdgeEvents), length(incorrectEdgeEvents)), " events)", sep=""))
      title(main=mainTitle)
      PV <- c()
      for (j in 1:7) {
        MF <- cbind(MC[,j], rowSums(MC[,-j]))
        fisherResult <- fisher.test(MF, alternative="two.sided")
        cat("Fisher test", eventNames[[j]], "vs branch correctness p-value =", fisherResult$p.value, "\n")
        PV <- c(PV, fisherResult$p.value)
      }
      text(x=colMeans(z), y=usrFromRelativeY(-0.1), xpd=TRUE, labels=paste0("p=", round(PV, 2)))
    }

    
    mainTitle <- "Event lengths (in number of genes)"
    correctEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect > 0]
    incorrectEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect < 0]
    rlcCompareHist(allEdgeEventsParams, ifelse(allEdgeEventsCorrect > 0, 1, 2),
                   relative=TRUE,
                   breaks=seq(-1, max(allEdgeEventsParams)), col=c("forestgreen","orange"),
                   add.legend=FALSE)
    legend("topright", fill=c("forestgreen", "orange"),
           legend=paste(c("correct splits", "incorrect splits"), " (", c(length(correctEdgeEvents), length(incorrectEdgeEvents)), " events)", sep=""))
    title(main=mainTitle)
    
    
    correctEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect > 0 & abs(allEdgeEventsCorrect) <= smallBranchThreshold]
    incorrectEdgeEvents <- allEdgeEvents[allEdgeEventsCorrect < 0 & abs(allEdgeEventsCorrect) <= smallBranchThreshold]
    keep <- abs(allEdgeEventsCorrect) <= smallBranchThreshold
    if (sum(keep) > 0) {
      rlcCompareHist(allEdgeEventsParams[keep], ifelse(allEdgeEventsCorrect[keep] > 0, 1, 2),
                     relative=TRUE,
                     breaks=seq(-1, max(allEdgeEventsParams)), col=c("forestgreen","orange"),
                     add.legend=FALSE)
      legend("topright", fill=c("forestgreen", "orange"),
             legend=paste(c("correct splits with <= ", "incorrect splits with <= "), smallBranchThreshold, " events", " (", c(length(correctEdgeEvents), length(incorrectEdgeEvents)), " events)", sep=""))
      title(main=mainTitle)
    }
    
  }
  
  invisible(dev.off())
  
  save.image(paste(params$outDirectory, "/benchmark-rdata-random.rda", sep=""))
  
  cat("Finished post-analysis.\n")
  
  sink()
}
