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
for (paramName in c("genes", "chr", "averageEvents", "nSpecies", "averageInversionLength", 
                    "averageDuplicationLength", "averageDeletionLength", "probInversion", 
                    "probTranslocation", "probDuplication", "probDeletion", "probFusion", 
                    "probFission", "probWGD", "averageDeletionRateWGD", "outDirectory", 
                    "nSimulations", "seed","minimalEventsOnBranch")) {
  if (is.null(params[[paramName]])) {
    stop(paste("Missing parameter:", paramName))
  }
}
rm(paramName)

#print(data.frame(Parameters=t(as.data.frame(params))))

if (!file.exists(params$outDirectory)) {
  dir.create(params$outDirectory)
}

# Check params
eventProbabilities <- c(params$probInversion, params$probTranslocation, params$probDuplication, params$probDeletion, params$probFusion, params$probFission, params$probWGD)
tmp <- sum(eventProbabilities)
if (tmp != 1) {
  cat("Probabilities sum to", tmp, " -> dividing by sum.\n")
  eventProbabilities <- eventProbabilities/sum(eventProbabilities)
}
rm(tmp)

# Make a PDF file that summarizes parameters
local({
  pdf(paste(params$outDirectory, "/parameters.pdf", sep=""))
  plot.new()
  title(main="Parameters")
  legend(0, 1, legend=c(names(params), as.character(params)), ncol=2, xpd=TRUE, box.lty=0)
  #x <- 1:50
  #y <- dpois(x, lambda=params$averageSpecies)
  #y[x < params$minimalSpecies] <- 0
  #plot(x, y, type="o", xlab="number of species", ylab="Probability", main="Number of species theoretical distribution")
  x <- 1:10
  plot(x, dpois(x, lambda=params$averageInversionLength), type="o", xlab="genes", ylab="Probability", main="Inversion length theoretical distribution")
  plot(x, dpois(x, lambda=params$averageDuplicationLength), type="o", xlab="genes", ylab="Probability", main="Duplication length theoretical distribution")
  plot(x, dpois(x, lambda=params$averageDeletionLength), type="o", xlab="genes", ylab="Probability", main="Deletion length theoretical distribution")
  dev.off()
})

# Run the actual simations
allSimOutDir <- paste(params$outDirectory, "/simulations", sep="")
if (!file.exists(allSimOutDir)) {
  dir.create(allSimOutDir)
}

if (!is.na(params$seed)) {
  set.seed(params$seed)
}

statsNSpecies <- rep(0, params$nSimulations)
statsNEvents <- matrix(0, nrow=params$nSimulations, ncol=7)
colnames(statsNEvents) <- eventNames
statsNEventsPerSpecies <- matrix(0, nrow=params$nSimulations*params$nSpecies, ncol=7)
colnames(statsNEventsPerSpecies) <- eventNames

#cat("Running ", params$nSimulations, " simulation", (if (params$nSimulations > 1) "s" else ""), "...\n", sep="")
#progressBarInit()

Rprof(NULL)
profFile <- paste(params$outDirectory, "/Rprof.out", sep="")
if (file.exists(profFile)) {
  invisible(file.remove(profFile))
}
#Rprof(profFile)
t0 <- proc.time()[["elapsed"]]

for (simIndex in 1:params$nSimulations) {
  cat("Running simulation ", simIndex, "/", params$nSimulations, " (t=", (proc.time()[["elapsed"]]-t0), ")...\n", sep="")
  
  simOutDir <- paste(allSimOutDir, "/", simIndex, sep="")
  if (!file.exists(simOutDir)) {
    dir.create(simOutDir)
  }
  
  # How many species are there?
  nSpecies <- as.integer(params$nSpecies)
  #while (nSpecies < params$minimalSpecies) {
  #  nSpecies <- rpois(1, lambda=params$averageSpecies)
  #}

  # When did the speciation events happen?
  speciationTimes <- sort(runif(n=nSpecies-1))
  # The root branch is useless, so assume the first speciation
  # is at the root.
  speciationTimes[[1]] <- 0
  
  # Generate tree
  tree <- 1
  for (s in seq(along=speciationTimes)) {
    tree <- splitRandomLeaf(tree, speciationTimes[[s]])
  }
  rm(s)
  
  # Convert tree to ape library format
  newickTree <- convertTreeToNewick(tree)
  apeTree <- read.tree(text=newickTree)
  apeTree$tip.label <- paste("sp", 1:length(apeTree$tip.label), sep="")
  rm(newickTree)
  
  # Generate the root genome
  genomes <- list()
  genomesStats <- list()
  rootIndex <- setdiff(apeTree$edge[1,], apeTree$edge[,2])
  stopifnot(length(rootIndex) == 1)
  local({
    genome <- list()
    genesChr <- tabulate(sample(1:params$chr, size=params$genes, replace=TRUE))
    # Remove empty chromosomes
    genesChr <- genesChr[genesChr != 0]
    stopifnot(sum(genesChr) == params$genes)
    k <- 1
    for (chr in 1:length(genesChr)) {
      genome[[chr]] <- addCentromere(seq.int(k, k+genesChr[[chr]]-1))
      k <- k + genesChr[[chr]]
    }
    genomes[[rootIndex]] <<- genome
  })
  genomesStats[[rootIndex]] <- list(eventStats=rep(0,7))
  
  # Simulate events on genomes
  events <- list()
  eventsParams <- list()
  eventsCount <- c()
  local({
    # Find correct order in which to simulate edges
    # (so that ancestors are computed before their children)
    orderedEdges <- findEdgeOrder(rootIndex)
    for (e in orderedEdges) {
      local({
        ancestor <- apeTree$edge[[e, 1]]
        child <- apeTree$edge[[e, 2]]
        genome <- genomes[[ancestor]]
        genomeStats <- genomesStats[[ancestor]]
        stopifnot(!is.null(genome))
        
        # How many events on branch?
        k <- as.integer(rpois(1, lambda=params$averageEvents*apeTree$edge.length[[e]]))
        if (k < params$minimalEventsOnBranch) {
          k<- as.integer(params$minimalEventsOnBranch)
        }
        eventsCount[[e]] <<- k
        
        edgeEvents <- rep(0L, k)
        edgeEventsParams <- rep(-1, k)
        for (evi in seq(along=edgeEvents)) {
          # What kind of event
          p <- eventProbabilities
          if (length(genome) < 2) {
            # Verify there is more that one chromosome,
            # otherwise some events are not possible.
            p[[2]] <- 0 # translocation
            p[[5]] <- 0 # fusion
          }
          ev <- sample.int(n=7L, size=1, replace=TRUE, prob=p)
          edgeEvents[[evi]] <- ev
          
          if (ev == 1L) {
            # Inversion
            chr <- randomChromosome(genome)
            l <- 0
            while (l < 1) { # retry if l == 0
              l <- rpois(1, lambda=params$averageInversionLength)
            }
            edgeEventsParams[[evi]] <- l
            r <- randomGeneSequence(genome[[chr]], l)
            genome[[chr]][r] <- reverseGenes(genome[[chr]][r])
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            stopifnot(all(!is.na(genome[[chr]][r])))
            rm(chr, l, r)
          } else if (ev == 2L) {
            # Reciprocal translocation
            chr1 <- randomChromosome(genome)
            chr2 <- randomChromosome(genome, exclude=chr1)
            # Find centromeres
            cm1 <- which(genome[[chr1]] == 0)
            cm2 <- which(genome[[chr2]] == 0)
            # Randomly selection positions on both chromosomes
            pos1 <- sample.int(n=length(genome[[chr1]])+1, size=1)
            pos2 <- sample.int(n=length(genome[[chr2]])+1, size=1)
            # Are we one the left or right of centromere?
            eventOnRight1 <- pos1 > cm1
            eventOnRight2 <- pos2 > cm2
            hasToReverseGenes <- eventOnRight1 != eventOnRight2
            s1 <- splitChromosome(genome[[chr1]], pos1)
            s2 <- splitChromosome(genome[[chr2]], pos2)
            if (eventOnRight1 && eventOnRight2) {
              newChr1 <- c(s1$left, s2$right)
              newChr2 <- c(s2$left, s1$right)
            } else if (eventOnRight1 && !eventOnRight2) {
              newChr1 <- c(s1$left, reverseGenes(s2$left))
              newChr2 <- c(reverseGenes(s1$right), s2$right)
            } else if (!eventOnRight1 && eventOnRight2) {
              newChr1 <- c(reverseGenes(s2$right), s1$right)
              newChr2 <- c(s2$left, reverseGenes(s1$left))
            } else { # both event are on the left of the centromere
              newChr1 <- c(s2$left, s1$right)
              newChr2 <- c(s1$left, s2$right)
            }
            # Check there is exactly one centromere on each new chromosome
            checkCentromere(newChr1)
            checkCentromere(newChr2)
            genome[[chr1]] <- newChr1
            genome[[chr2]] <- newChr2
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            rm(chr1, chr2, cm1, cm2, pos1, pos2, eventOnRight1, eventOnRight2,
               hasToReverseGenes, s1, s2, newChr1, newChr2)
          } else if (ev == 3L) {
            # Duplication
            chr <- randomChromosome(genome)
            l <- 0
            while (l < 1) { # retry if l == 0
              l <- rpois(1, lambda=params$averageDuplicationLength)
            }
            edgeEventsParams[[evi]] <- l
            r <- randomGeneSequence(genome[[chr]], l)
            duplicatedSequence1 <- genome[[chr]][r]
            duplicatedSequence2 <- duplicatedSequence1
            if (any(duplicatedSequence1 == 0)) {
              # We have the centromere, take care not to duplicate it.
              # Select at random one of the two copies to have the centromere.
              if (runif(1) > 0.5) {
                duplicatedSequence1 <- duplicatedSequence1[duplicatedSequence1 != 0]
              } else {
                duplicatedSequence2 <- duplicatedSequence2[duplicatedSequence2 != 0]
              }
            }
            genome[[chr]] <- c(genome[[chr]][a2b(1, r[[1]]-1)], 
                               duplicatedSequence1,
                               duplicatedSequence2,
                               genome[[chr]][a2b(r[[length(r)]]+1, length(genome[[chr]]))])
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            # Check there is exactly one centromere
            checkCentromere(genome[[chr]])
            rm(chr, l, r, duplicatedSequence1, duplicatedSequence2)
          } else if (ev == 4L) {
            # Deletion
            chr <- randomChromosome(genome)
            l <- 0
            while (l < 1) { # retry if l == 0
              l <- rpois(1, lambda=params$averageDeletionLength)
            }
            edgeEventsParams[[evi]] <- l
            r <- randomGeneSequence(genome[[chr]], l)
            # Do not remove the centromere
            r <- r[genome[[chr]][r] != 0]
            genome[[chr]] <- genome[[chr]][-r]
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            # Check there is exactly one centromere
            checkCentromere(genome[[chr]])
            rm(chr, l, r)
          } else if (ev == 5L) {
            # Fusion
            chr1 <- randomChromosome(genome)
            chr2 <- randomChromosome(genome, exclude=chr1)
            genesChr1 <- genome[[chr1]]
            genesChr2 <- genome[[chr2]]
            # Relative orientation of the 2 chromosomes to be fused
            if (runif(1) > 0.5) {
              genesChr2 <- reverseGenes(genesChr2)
            }
            if (runif(1) > 0.5) {
              newChr <- c(genesChr1, removeCentromere(genesChr2))
            } else{
              newChr <- c(removeCentromere(genesChr1), genesChr2)
            }
            stopifnot(sum(newChr== 0) == 1)
            genome[[chr1]] <- newChr
            genome[[chr2]] <- NULL
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            rm(chr1, chr2, genesChr1, genesChr2, newChr)
          } else if (ev == 6L) {
            # Fission
            chr <- randomChromosome(genome)
            pos <- sample.int(n=length(genome[[chr]])+1, size=1)
            s <- splitChromosome(genome[[chr]], pos)
            if (sum(s$left == 0) == 0) {
              s$left <- addCentromere(s$left)
            } else {
              s$right <- addCentromere(s$right)
            }
            genome[[chr]] <- s$left
            genome[[length(genome)+1]] <- s$right
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            rm(chr, pos, s)
          } else if (ev == 7L) {
            # Whole Genome Duplication
            nchr <- length(genome)
            for (chr in 1:nchr) {
              cm <- findCentromere(genome[[chr]])
              newChr1 <- genome[[chr]]
              newChr2 <- genome[[chr]]
              deleteOneCopy <- runif(n=length(genome[[chr]])) < params$averageDeletionRateWGD
              deleteOneCopy[[cm]] <- FALSE # never delete the centromere
              chooseCopy1 <- runif(sum(deleteOneCopy)) > 0.5
              deleteCopy1 <- deleteOneCopy
              deleteCopy2 <- deleteOneCopy
              deleteCopy1[deleteOneCopy] <- chooseCopy1
              deleteCopy2[deleteOneCopy] <- !chooseCopy1
              newChr1 <- newChr1[!deleteCopy1]
              newChr2 <- newChr2[!deleteCopy2]
              checkCentromere(newChr1)
              checkCentromere(newChr2)
              # Remove a chromosome if it contains no gene
              if (length(newChr1) == 1) newChr1 <- NULL
              if (length(newChr2) == 1) newChr2 <- NULL
              genome[[chr]] <- newChr1
              genome[[length(genome)+1]] <- newChr2
            }
            genomeStats$eventStats[[ev]] <- genomeStats$eventStats[[ev]] + 1
            rm(cm, newChr1, newChr2, deleteOneCopy, chooseCopy1, deleteCopy1, deleteCopy2)
          } else {
            stop("bug")
          }
          nGenesOnChr <- sapply(genome, length) - 1L
          if (any(nGenesOnChr == 0)) {
            genome <- genome[-which(nGenesOnChr == 0)]
          }
        }
        
        events[[e]] <<- edgeEvents
        eventsParams[[e]] <<- edgeEventsParams
        genomes[[child]] <<- genome
        genomesStats[[child]] <<- genomeStats
      })
    }
  })
  
  
  
  # Write the genomes in GRIMM-like format.
  # (Not exactly GRIMM format because GRIMM cannot handle duplications and deletions.)
  local({
    fh <- file(paste(simOutDir, "/genomes.txt", sep=""), open="w")
    writeLines(">ancestor", con=fh)
    writeGenome(genomes[[rootIndex]], fh)
    for (species in 1:nSpecies) {
      genome <- genomes[[species]]
      writeLines(paste(">", apeTree$tip.label[[species]], sep=""), con=fh)
      writeGenome(genomes[[species]], fh)
    }
    close(fh)
  })
  
  
  # Write tree in Newick format (branch length = time)
  write.tree(apeTree, file=paste(simOutDir, "/tree-newick-time.txt", sep=""))
  
  # Write tree in Newick format (branch length = number of events)
  apeTreeE <- apeTree
  apeTreeE$edge.length <- eventsCount
  write.tree(apeTreeE, file=paste(simOutDir, "/tree-newick.txt", sep=""))
  
  # Draw the tree
  local({
    pdf(paste(simOutDir, "/tree.pdf", sep=""), width=11, height=8)
    
    generalTitle <- paste("\nsimulation ", simIndex, "/", params$nSimulations, sep="")
    
    chrBySpecies <- unlist(lapply(genomes[1:nSpecies], length))
    genesBySpecies <- rep(0, nSpecies)
    for (s in 1:nSpecies) {
      nGenesOnChr <- sapply(genomes[[s]], length) - 1L
      genesBySpecies[[s]] <- sum(nGenesOnChr)
    }
    apeTree$tip.label <- paste(apeTree$tip.label, " c=", chrBySpecies, " g=", genesBySpecies, sep="")
    pl <- plot.phylo(apeTree)
    title(main=paste("Events on each branch", generalTitle, sep=""))
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    for (e in 1:nrow(apeTree$edge)) {
      if (eventsCount[[e]] > 0) {
        ancestor <- apeTree$edge[[e,1]]
        child <- apeTree$edge[[e,2]]
        x0 <- lastPP$xx[ancestor]
        y0 <- lastPP$yy[child]
        x1 <- lastPP$xx[child]
        y1 <- lastPP$yy[child]
        x <- x0 + (x1-x0)*seq(0, 1, length.out=eventsCount[[e]]+2)[2:(eventsCount[[e]]+1)]
        y <- rep(y0, eventsCount[[e]])
        points(x, y, col=eventColors[events[[e]]], pch=eventPch[events[[e]]], cex=eventCex[events[[e]]])
      }
    }
    legend(x=usrFromRelativeX(-0.05), y=usrFromRelativeY(-0.05), xpd=TRUE, legend=eventNames, pch=eventPch, col=eventColors, horiz=TRUE, pt.cex=eventCex)
    
    pl <- plot.phylo(apeTree)
    title(main=paste("Total number of events on each branch", generalTitle, sep=""))
    edgelabels(eventsCount)
    
    plot.phylo(apeTreeE)
    title(main=paste("Total number of events on each branch (x axis = events)", generalTitle, sep=""))
    edgelabels(eventsCount)
    
    dev.off()
  })
  
  
  
  # Write in PhyChro format
  if (!is.null(params$phyChroFormat) && params$phyChroFormat) {
    phyChroDir <- paste(simOutDir, "/PhyChro", sep="")
    if (!file.exists(phyChroDir)) {
      dir.create(phyChroDir)
    } else {
      # Clear directory contents
      invisible(file.remove(list.files("results_yeast/simulations/3/PhyChro", full.names=TRUE, pattern="^S")))
    }
    if (nSpecies > 999) {
      stop("too many (>999) species for PhyChro format")
    }
    speciesPhyChroNames <- sprintf("S%03d", 1:nSpecies)
    phyChroDefAll <- list()
    for (species in 1:nSpecies) {
      local({
        genome <- genomes[[species]]
        speciesPhyChroName <- speciesPhyChroNames[[species]]
        featuresByChr <- sapply(genome, length)
        nFeatures <- sum(featuresByChr)
        featureChr <- as.integer(rep(1:length(genome), featuresByChr))
        stopifnot(length(featureChr) == nFeatures)
        allFeatures <- unlist(genome)
        stopifnot(length(allFeatures) == nFeatures)
        IDg.chr <- rep(0L, nFeatures)
        IDg.all <- rep(0L, nFeatures)
        IDf.all <- rep(0L, nFeatures)
        phyChroCh <- matrix(0L, nrow=3, ncol=length(genome))
        phyChroCh[1,] <- 1:length(genome)
        phyChroCh[2,] <- which(allFeatures == 0)
        for (chr in 1:length(genome)) {
          IDg.chr[featureChr == chr] <- cumsum(genome[[chr]] != 0)
          phyChroCh[[3,chr]] <- as.integer(max(which(featureChr == chr)))
        }
        IDg.all <- cumsum(allFeatures != 0)
        phyChroDef <- data.frame(type=ifelse(allFeatures != 0, "gene", "centromere"),
                                 name=paste("sp", species, "f", as.integer(1:nFeatures), sep=""),
                                 chr=as.integer(featureChr),
                                 start=rep(0L, nFeatures),
                                 end=rep(0L, nFeatures),
                                 strand=rep(0L, nFeatures),
                                 sens=rep(0L, nFeatures),
                                 "IDg/chr"=as.integer(IDg.chr),
                                 "IDg/all"=as.integer(IDg.all),
                                 "IDf/all"=as.integer(1:nFeatures),
                                 check.names=FALSE)
        phyChroDef$type[allFeatures == 0] <- "centromere"
        
        write.table(phyChroCh, file=paste(phyChroDir, "/", speciesPhyChroName, ".ch", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        write.table(phyChroDef, file=paste(phyChroDir, "/", speciesPhyChroName, ".def", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
        
        phyChroDefAll[[species]] <<- phyChroDef
      })
    }
    rm(species)
    
    for (speciesA in a2b(1L, nSpecies-1L)) {
      genomeA <- genomes[[speciesA]]
      genomeAflat <- unlist(genomeA)
      A <- phyChroDefAll[[speciesA]]
      A$homolog <- abs(genomeAflat)
      A$forward <- genomeAflat > 0
      
      for (speciesB in a2b(speciesA+1L, nSpecies)) {
        genomeB <- genomes[[speciesB]]
        genomeBflat <- unlist(genomeB)
        B <- phyChroDefAll[[speciesB]]
        B$homolog <- abs(genomeBflat)
        B$forward <- genomeBflat > 0
        
        # Make a matrix of all homologous pairs
        M <- merge(A, B, by="homolog")
        colnames(M) <- gsub("/", ".", colnames(M))
        # Remove centromere
        M <- M[M$homolog != 0, ]
        # Order matrix by order of first genome
        M <- M[order(M$IDf.all.x), ]
        M$sameStrand <- M$forward.x == M$forward.y
        
        # Compute synteny blocks
        # The super-slow part of the program
        #system.time({
        M$syntenyBlock <- 1:nrow(M)
        for (chrA in 1:length(genomeA)) {
          for (chrB in 1:length(genomeB)) {
            rowIndices <- which(M$chr.x == chrA & M$chr.y == chrB)
            if (length(rowIndices) != 0) {
              subM <- M[rowIndices, c("IDg.all.x","IDg.all.y","syntenyBlock","sameStrand")]
              for (i in 1:nrow(subM)) {
                if (subM$sameStrand[[i]]) {
                  subM$syntenyBlock[subM$IDg.all.x == subM$IDg.all.x[[i]] + 1L & subM$IDg.all.y == subM$IDg.all.y[[i]] + 1L] <- subM$syntenyBlock[[i]]
                } else {
                  subM$syntenyBlock[subM$IDg.all.x == subM$IDg.all.x[[i]] + 1L & subM$IDg.all.y == subM$IDg.all.y[[i]] - 1L] <- subM$syntenyBlock[[i]]
                }
              }
              M$syntenyBlock[rowIndices] <- subM$syntenyBlock
            }
          }
        }
        #})
        
        # Remove blocks with only one pair
        blocksOK <- duplicated(M$syntenyBlock) | duplicated(M$syntenyBlock, fromLast=TRUE)
        M <- M[blocksOK, ]

        # Write synteny blocks
        synBlocks <- unique(M$syntenyBlock)
        sbLines <- rep(NA_character_, length(synBlocks))
        if (length(synBlocks) != 0) {
          for (i in 1:length(synBlocks)) {
            w <- which(M$syntenyBlock == synBlocks[[i]])
            sbString <- paste(M$chr.x[[w[[1]]]], M$chr.y[[w[[1]]]])
            for (g in w) {
              sbString <- paste(sbString, paste(M[g,c("name.x","IDf.all.x","name.y","IDf.all.y")], collapse=" "), 100L)
            }
            sbLines[[i]] <- sbString
          }
        }
        writeLines(sbLines, con=paste(phyChroDir, "/", speciesPhyChroNames[[speciesA]], ".", speciesPhyChroNames[[speciesB]], ".orth.synt", sep=""))
        write.table(data.frame(name1=M$name.x,
                               chr1=M$chr.x,
                               IDr.all1=rep(0L, nrow(M)),
                               IDg.all1=M$IDg.all.x,
                               IDf.all1=M$IDf.all.x,
                               name2=M$name.y,
                               chr2=M$chr.y,
                               IDr.all2=rep(0L, nrow(M)),
                               IDg.all2=M$IDg.all.y,
                               IDf.all2=M$IDf.all.y,
                               simi=rep(100L, nrow(M)),
                               sstr=ifelse(M$sameStrand, "sameW", "invert")),
                    file=paste(phyChroDir, "/", speciesPhyChroNames[[speciesA]], ".", speciesPhyChroNames[[speciesB]], ".orth.pairs", sep=""), quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
        

      
      } # end species B
    } # end species A
  } # end PhyChro format
  
  # Save all data
  if (!is.null(params$save) && params$save) {
    save(simIndex, genomes, genomesStats, apeTree, apeTreeE, events, eventsCount, eventsParams,
         file=paste(simOutDir, "/sim-R-data.RData", sep=""))
  }
  
  # Record stats
  statsNSpecies[[simIndex]] <- nSpecies
  statsNEvents[simIndex,] <- tabulate(unlist(events), nbins=7)
  
  statsNEventsPerSpecies[ params$nSpecies*(simIndex-1)+seq(1,nSpecies) ,] <- t(sapply(genomesStats[1:nSpecies], function(gs) unlist(gs$eventStats)))
  
  # Show progress bar
  #progressBarNext(simIndex/params$nSimulations)
}

t1 <- proc.time()[["elapsed"]]
Rprof(NULL)
#progressBarEnd()
simTime <- t1 - t0
cat("Total simulation time: ", simTime, " s (", (simTime/params$nSimulations), " s per simulation)\n", sep="")

pdf(paste(params$outDirectory, "/statistics.pdf", sep=""))
if (params$nSimulations > 1) {
  invisible(rlcCompareDensMatrix(statsNEvents, main="Total number of events in trees", col=eventColors))
  #plot(density(statsNSpecies), main=paste("Number of species\nmean =", round(mean(statsNSpecies), 1), " param =", params$averageSpecies))
}
for (i in 1:length(eventNames)) {
  hist(statsNEvents[,i], main=paste(eventNames[[i]], " (mean = ", mean(statsNEvents[,i]), ")", sep=""))
}
invisible(dev.off())


pdf(paste(params$outDirectory, "/statistics per species.pdf", sep=""))
tmp <- rowSums(statsNEventsPerSpecies)
invisible(hist(tmp, main=paste("Total number of events from ancestor to species\nmean = ", mean(tmp), sep=""), xlab="Number of events"))
if (nrow(statsNEventsPerSpecies) > 1) {
  invisible(rlcCompareDensMatrix(statsNEventsPerSpecies, main="Number of events from ancestor to species", col=eventColors))
}
for (i in 1:length(eventNames)) {
  hist(statsNEventsPerSpecies[,i], main=paste(eventNames[[i]], " (mean = ", mean(statsNEventsPerSpecies[,i]), ")", sep=""))
}
invisible(dev.off())


cat("Finished.\n")

# Print profiling information
if (file.exists(profFile)) {
  sink(paste(params$outDirectory, "/Rprof.txt", sep=""))
  print(summaryRprof(paste(params$outDirectory, "/Rprof.out", sep="")))
  sink(NULL)
}


