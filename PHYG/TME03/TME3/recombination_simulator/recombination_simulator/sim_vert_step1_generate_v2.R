# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

setwd("~/Documents/PhyChro/sim")

params <- list(
  genes = 18000,
  chr = 23,
  nSpecies = 13,
  averageEvents = 1000,
  minimalEventsOnBranch = 1,
  averageInversionLength = 5,
  averageDuplicationLength = 5,
  averageDeletionLength = 1,
  probInversion = 0.6,
  probTranslocation = 0.2979,
  probDuplication = 0.05,
  probDeletion = 0.05,
  probFusion = 0.001,
  probFission = 0.001,
  probWGD = 0.0001,
  averageDeletionRateWGD = 0.8,
  outDirectory = "results_vert1",
  #phyChroFormat = TRUE,
  save = TRUE,
  nSimulations = 100,
  seed = 0)

source("scripts/sim_main.R")
