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
for (paramName in c("outDirectory","phyChroInstallationDirectory")) {
  if (is.null(params[[paramName]])) {
    stop(paste("Missing parameter:", paramName))
  }
}
rm(paramName)

# Translate to absolute path: safer because PhyChro requires we change working directory
originalWD <- getwd()
allSimOutDir <- normalizePath(paste(params$outDirectory, "/simulations", sep=""))

simIndex <- 1L
t0 <- proc.time()[["elapsed"]]
while (TRUE) {
  simOutDir <- paste(allSimOutDir, "/", simIndex, sep="")
  phyChroDir <- paste(simOutDir, "/PhyChro", sep="")
  if (!file.exists(phyChroDir)) {
    break
  }
  
  cat("Running PhyChro for simulation", simIndex, "...\n")
  
  setwd(params$phyChroInstallationDirectory)
  ret <- system2("python", args=c("PhyChro.py", paste(phyChroDir, "/", sep=""), "tree", "0"))
  setwd(originalWD)
  if (ret != 0) {
    stop("PhyChro failed")
  }
  
  phyChroTree <- read.tree(paste(phyChroDir, "/tree.outtree", sep=""))
  correctTree <- read.tree(paste(simOutDir, "/tree-newick.txt", sep=""))
  correctTreeUnrooted <- unroot(correctTree)
  
  # Translate back species names to spXXX
  phyChroTree$tip.label <- paste("sp", as.integer(gsub("S", "", phyChroTree$tip.label)), sep="")
  
  pdf(paste(simOutDir, "/PhyChroTree.pdf", sep=""), width=11, height=8)
  
  layout(matrix(1:2, ncol=2))
  
  plot.phylo(phyChroTree, type="unrooted")
  title(main="PhyChro tree")
  plot.phylo(correctTreeUnrooted, type="unrooted")
  title(main="Correct tree")
  
  plot.phylo(phyChroTree, use.edge.length=FALSE, type="unrooted")
  title(main="PhyChro tree (no length)")
  plot.phylo(correctTreeUnrooted, use.edge.length=FALSE, type="unrooted")
  title(main="Correct tree (no length)")
  
  dev.off()
  
  simIndex <- simIndex + 1L
}

setwd(originalWD)

t1 <- proc.time()[["elapsed"]]
simTime <- t1 - t0
cat("Total simulation time: ", simTime, " s (", (simTime/(simIndex-1)), " s per simulation)\n", sep="")



  