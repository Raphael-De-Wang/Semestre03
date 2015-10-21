# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

options(stringsAsFactors=FALSE)

source("scripts/auxfunctions.R")

library(ape)

# inv trans dup del fus fiss WGD
eventNames <- c("inversion", "translocation", "duplication", "deletion", "fusion", "fission", "WGD")
eventShortNames <- c("inv", "trans", "dup", "del", "fus", "fis", "WGD")
eventColors <- c("black", "blue", "green", "orange", "magenta", "darkgreen", "red")
eventPch <- c(39, 4, 2, 6, 7, 8, 5)
eventCex <- c(0.8, 0.8, 1, 1, 2, 2, 3)

splitRandomLeaf <- function(node, t) {
  if (length(node) > 1) {
    chosenChildIndex <- sample.int(2, size=1)
    child <- node[[chosenChildIndex]]
    return(c(list(splitRandomLeaf(child, t)), node[-chosenChildIndex]))
  } else {
    l <- node
    stopifnot(is.numeric(l))
    ls <- (1-t)
    stopifnot(ls >= 0)
    la <- l - ls
    ancestor <- list(ls, ls)
    ancestor$length <- la
    return(ancestor)
  }
}

convertTreeToNewickAux <- function(tree) {
  if (is.list(tree)) {
    stopifnot(length(tree) > 1)
    return(paste("(", convertTreeToNewickAux(tree[[1]]), ",", convertTreeToNewickAux(tree[[2]]), "):", tree$length, sep=""))
  } else {
    return(paste("species:", tree, sep=""))
  }
}


convertTreeToNewick <- function(tree) {
  return(paste(convertTreeToNewickAux(tree), ";", sep=""))
}


# Functions to display progress bars
progressBarInit <- function() {
  for (i in 1:10) {
    cat(" ")
  }
  for (i in 1:50) {
    cat("_")
  }
  cat("\n")
  for (i in 1:10) {
    cat(" ")
  }
  .progressBarLast <<- 0
}

# x is in [0,1]
progressBarNext <- function(x) {
  if (length(x) == 0) {
    stop("progressBarNext: invalid argument")
  }
  if (x > 1) x <- 1
  while (50*x >= .progressBarLast + 1) {
    cat("X")
    .progressBarLast <<- .progressBarLast + 1
  }
}

progressBarEnd <- function() {
  while (.progressBarLast < 50) {
    cat("X")
    .progressBarLast <<- .progressBarLast + 1
  }
  cat(" done.\n")
}

# For the current plot, takes a position as a number in [0,1]
# and return the position as units of the current plot.
# Example: if range is [100,200], passing 0.5 gives 150.
# Also work for <0 and >1 for out-of-plot things.
usrFromRelativeX <- function(x) {
  usr <- par("usr")
  return(x * (usr[2] - usr[1]) + usr[1])
}

# Same for Y
usrFromRelativeY <- function(y) {
  usr <- par("usr")
  return(y * (usr[4] - usr[3]) + usr[3])
}


findEdgeOrder <- function(nodeIndex) {
  edgesToChildren <- which(apeTree$edge[,1] == nodeIndex)
  if (length(edgesToChildren) > 0) {
    stopifnot(length(edgesToChildren) == 2)
    return(c(edgesToChildren, findEdgeOrder(apeTree$edge[[edgesToChildren[[1]],2]]), findEdgeOrder(apeTree$edge[[edgesToChildren[[2]],2]])))
  } else {
    return(NULL)
  }
}

# Select a chromosome randomly with a probability proportional
# to its number of genes.
# If exclude in specified (integer), it means the corresponding
# chromosome cannot be chosen.
randomChromosome <- function(genome, exclude=NULL) {
  chrProb <- sapply(genome, length) - 1 # -1 to remove centromere
  if (!is.null(exclude)) {
    chrProb[[exclude]] <- 0
  }
  return(sample.int(n=length(genome), size=1, prob=chrProb))
}

# Select a random sequence of k consecutive genes.
# If chromosome has less than k genes, k is automatically
# reduced to the number of genes.
# Note that the returned list my be of length k or k+1,
# because it can contain the centromere
# (which does not count as a gene).
# The function returns the indices of the elements in chromosome,
# ie. if we call r the returned value,
# access the genes using chromosome[r].
randomGeneSequence <- function(chromosome, k) {
  stopifnot(k > 0)
  if (k > length(chromosome) - 1) {
    k <- length(chromosome) - 1
  }
  k <- as.integer(k)
  a <- sample.int(n=length(chromosome)-k, size=1L)
  b <- a + k - 1L
  if (sum(chromosome[a:b] == 0) == 1) {
    # We have included the centromere, so add an extra gene to get the right count
    b <- b + 1
  }
  r <- seq(a,b)
  stopifnot(sum(chromosome[r] != 0) == k)
  return(r)
}

# Range from integer a to integer b.
# Example: a2b(4, 7) returns c(4L, 5L, 6L, 7L)
# It differs from seq.int(a,b) by the fact that
# if a>b, it returns an empty vector.
a2b <- function(a, b) {
  if (a <= b) {
    seq.int(a, b)
  } else {
    integer(0)
  }
}

# Split chromosome at given position
# Example: c(1, 2, 3, 4, 5, 6)  at position 4 will be:
# left=c(1, 2, 3)  right=c(4,5,6)
splitChromosome <- function(chromosome, position) {
  stopifnot(position >= 1 && position <= length(chromosome)+1)
  list(left=chromosome[a2b(1, position-1)], right=chromosome[a2b(position, length(chromosome))])
}

# Reverse a piece of DNA
reverseGenes <- function(genomeSequence) {
  # If a centromere is present (0), it will be left as is (since -0 == 0)
  - rev(genomeSequence)
}

checkCentromere <- function(chromosome) {
  cm <- which(chromosome == 0)
  if (length(cm) > 1) stop("more than one centromere")
  if (length(cm) == 0) stop("missing centromere")
  invisible(NULL)
}

findCentromere <- function(chromosome) {
  cm <- which(chromosome == 0)
  if (length(cm) > 1) stop("more than one centromere")
  if (length(cm) == 0) stop("missing centromere")
  return(cm)
}

# Remove centromere from chromosome
removeCentromere <- function(chromosome) {
  chromosome[chromosome != 0]
}

# Add centromere randomly on chromosome
addCentromere <- function(chromosome) {
  # Select a position between two genes to be the centromere
  centromere <- sample.int(n=length(chromosome)+1, size=1)
  # Add centromere (marked by 0)
  tmp <- splitChromosome(chromosome, centromere)
  return(as.integer(c(tmp$left, 0L, tmp$right)))
}


writeGenome <- function(genome, fh) {
  for (chr in 1:length(genome)) {
    chrGenome <- genome[[chr]]
    chrGenome <- chrGenome[chrGenome != 0]
    chrAsString <- paste(paste(chrGenome, collapse=" "), "$")
    writeLines(chrAsString, con=fh)
  }
}

eventsAsChar <- function(events, eventsParams) {
  r <- eventShortNames[events]
  s <- eventsParams >= 0
  r[s] <- paste(r[s], "(", eventsParams[s], ")", sep="")
  return(r)
}
