# Author: Raphael Champeimont
# UMR 7238 Biologie Computationnelle et Quantitative

# setwd("~/Documents/PhyChro/sim")

for (resName in list.files(pattern="^results")) {
  if (file.exists(paste(resName, "/simulations", sep=""))) {
    A <- c("sim_any_step3_analyze.R", shQuote(resName))
    cat(paste(A, collapse=" "), "\n")
    stopifnot(system2("Rscript", args=A) == 0)
    a <- paste(resName, "/plot-benchmark.pdf", sep="")
    if (file.exists(a)) {
      b <- paste("summary/", resName, ".pdf", sep="")
      stopifnot(file.copy(a, b, overwrite=TRUE))
    }
  }
}
