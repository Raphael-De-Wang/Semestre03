# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

setwd("~/Documents/PhyChro/sim")

params <- list(
  genes = 500,
  chr = 4,
  nSpecies = 5,
  averageEvents = 50,
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
  outDirectory = "results_test",
  phyChroFormat = TRUE,
  nSimulations = 5,
  seed = 0,
  save = TRUE,
  phyChroInstallationDirectory = "~/Documents/Software/CHROnicle/Programs/2PhyChro")

source("scripts/sim_main.R")
source("scripts/sim_run_phychro.R")
source("scripts/sim_analyze_phychro.R")
