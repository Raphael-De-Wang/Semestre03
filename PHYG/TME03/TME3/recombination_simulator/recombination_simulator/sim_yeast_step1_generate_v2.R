# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

# setwd("~/Documents/PhyChro/sim")

params <- list(
  genes = 5000,
  chr = 8,
  nSpecies = 10,
  averageEvents = 400,
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
  outDirectory = "results_yeast1",
  phyChroFormat = TRUE,
  phyChroInstallationDirectory = "/users/nfs/Etu9/3404759/Workspace/Semestre03/PHYG/TME03/2PhyChro/",
  save = TRUE,
  nSimulations = 10,
  seed = 0)

source("scripts/sim_main.R")
