# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

outdir <- commandArgs(trailingOnly=TRUE)[[1]]

setwd("~/Documents/PhyChro/sim")

stopifnot(file.exists(outdir))

source("~/R-scripts/auxfunctions.R")
source("~/R-scripts/auxlgm.R")

for (simIndex in 1:100) {
  script <- "cd ~/Documents/PhyChro/sim"
  script <- paste(script, paste("Rscript", "sim_any_step2_phychro_selected.R", outdir, simIndex), sep="\n")
  runOnLGMcluster(jobName=paste0("sim", simIndex), scriptContents=script)
}
