# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

setwd("~/Documents/PhyChro/sim")

params <- list(
  genes = 5000,
  chr = 8,
  nSpecies = 21,
  averageEvents = 500,
  minimalEventsOnBranch = 10,
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
  outDirectory = "results_yeast10",
  #phyChroFormat = TRUE,
  save = TRUE,
  nSimulations = 100,
  seed = 0)

source("scripts/sim_main.R")
