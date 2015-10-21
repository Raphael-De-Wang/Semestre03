# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

outdir <- commandArgs(trailingOnly=TRUE)[[1]]
simIndex <- as.integer(commandArgs(trailingOnly=TRUE)[[2]])

setwd("~/Documents/PhyChro/sim")

params <- list(
  outDirectory = outdir,
  phyChroInstallationDirectory = "~/Documents/Software/CHROnicle/Programs/2PhyChro",
  simIndex = simIndex)

source("scripts/sim_run_phychro_selected.R")
