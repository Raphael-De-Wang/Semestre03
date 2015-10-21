# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

outdir <- commandArgs(trailingOnly=TRUE)[[1]]

setwd("~/Documents/PhyChro/sim")

params <- list(
  outDirectory = outdir)

source("scripts/sim_analyze_random.R")
