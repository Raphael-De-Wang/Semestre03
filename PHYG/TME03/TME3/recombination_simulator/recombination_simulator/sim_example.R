# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

#setwd("~/Documents/PhyChro/sim")

params <- list(
  # Genes in ancestor:
  genes = 100,
  
  # Chromosomes in ancestor:
  chr = 4,
  
  # Numbef of species
  nSpecies = 5,
  
  # Average number of events from ancestor to final species:
  averageEvents = 5,
  
  # Normally the number of events on a branch is computed as the result
  # of a Poisson distribution of parameter lambda = branch length x averageEvents.
  # However, this may result in the number of events to be very small if
  # two speciation events are close to each other. If it is less than
  # minimalEventsOnBranch, the number of events is replaced by minimalEventsOnBranch.
  # This is usefull to benchmark tree inference methods because no events
  # or very few might make the task "impossible" for the program.
  # Note: This will make the average number of events a bit higher than averageEvents.
  # Set to 0 to disable this feature.
  minimalEventsOnBranch = 1,
  
  # Several length parameters (unit = number of genes)
  # Note that a value has no effect if the corresponding
  # event probability (see below) is 0.
  averageInversionLength = 5,
  averageDuplicationLength = 5,
  averageDeletionLength = 1,
  
  # Probability of each kind of event:
  # It is OK if they do not sum to 1, they will automatically
  # be divided by their sum.
  probInversion = 0.6,
  probTranslocation = 0.25,
  probDuplication = 0,
  probDeletion = 0,
  probFusion = 0.05,
  probFission = 0.05,
  probWGD = 0.05,
  
  # If a WGD occurs, for each gene, either one copy is deleted, or both are retained.
  # This is the probability that the first case occurs.
  # If set to 0, all copies are retained (genome size is doubled).
  # If set to 1, the number of copies of each gene remains identical
  # (genome size is constant).
  # This parameters has has no effect if probWGD is 0.
  averageDeletionRateWGD = 1,
  
  # Output directory:
  outDirectory = "results_example",
  
  # Whether or not to create additionnal files in PhyChro format.
  # (PhyChro: http://www.lcqb.upmc.fr/CHROnicle/)
  # Note that this will take much longer (like 10 times) if enabled,
  # because PhhyChro requires we compute synteny blocks.
  # If set to FALSE, only the genomes.txt and tree.pdf files are created.
  phyChroFormat = TRUE,
  
  # Number of simulations to make
  nSimulations = 10,
  
  # Random seed:
  # Set to NA to produce a different set of simulations every time your run the script.
  # Set to a integer to be able to reproduce the SAME set of simulations again
  # (provided of course you don't change the parameters).
  seed = 0,
  
  # Save R data in a file for post-analysis.
  save = TRUE
  )

source("scripts/sim_main.R")
