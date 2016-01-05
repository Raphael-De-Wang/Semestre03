#!env Rscript

# 1. Download the file MetagenomicsExp.R containing
temp.space <- new.env()
load("MetagenomicsExp.R",temp.space)
# ls(temp.space)
gene_ab <- temp.space$gene_ab
y <- temp.space$y

# 2. Install the kernlab R package
library(kernlab)

# 3. tutorial on SVMs with kernlab
# http://cbio.ensmp.fr/~jvert/svn/tutorials/practical/svmbasic/svmbasic_notes.pdf
# source("svm_intro.R")

# 4. Train a linear SVM on the training (metagenomics) set
svp <- ksvm(gene_ab,y,type="C-svc",kernel='vanilladot',C=100,scaled=c())
