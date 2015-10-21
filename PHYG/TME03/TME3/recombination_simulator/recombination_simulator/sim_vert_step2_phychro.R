# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

setwd("~/Documents/PhyChro/sim")

params <- list(
  outDirectory = "results_vert",
  phyChroInstallationDirectory = "~/Documents/Software/CHROnicle/Programs/2PhyChro")

source("scripts/sim_run_phychro.R")
