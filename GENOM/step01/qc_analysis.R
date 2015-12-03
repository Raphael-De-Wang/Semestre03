#!/usr/bin/Rscript

require(optparse)

option_list = list(
  make_option(c("--lmiss"), action="store", default=NULL, type='character',
              help="MISSING RATES"),
  make_option(c("--imiss"), action="store", default=NULL, type='character',
              help="MISSING RATES"),
  make_option(c("--hardy"), action="store", default=NULL, type='character',
              help="Hardy-Weinberg Equilibrium"),
  make_option(c("--maf"),   action="store", default=NULL, type='character',
              help="MAF - Minor allele frequency"),
  make_option(c("--save"),  action="store", default=NULL, type='character',
              help="save plot")  
)


opt = parse_args(OptionParser(option_list=option_list))


lmiss_call_rate <- function(fname, taille=50) {
    lmiss <- list(data=NA,x=NA,y=NA)
    lmiss$data <- read.table(fname,header=T)
    lmiss$x <- seq(0.91,0.99, length.out=taille)
    lmiss$y <- c()
    for (i in 1:taille){
	lmiss$y <- c(sum(lmiss$data[,5]<=(1-lmiss$x[i]))/length(lmiss$data[,5]),lmiss$y)
    }
    plot(lmiss$y~lmiss$x, type='l', col='blue', main="SNP")
    return(lmiss)
}


imiss_call_rate <- function(fname, taille=50) {
    imiss <- list(data=NA,x=NA,y=NA)
    imiss$data <- read.table(fname, header=T)
    imiss$x <- seq(0.91,0.99, length.out=taille)
    imiss$y <- c()
    for (i in 1:taille){
	imiss$y <- c(sum(imiss$data[,6]<=(1-imiss$x[i]))/length(imiss$data[,6]),imiss$y)
    }
    plot(imiss$y~imiss$x, type='l', col='red', main="FID")
    return(imiss)
}


hardy_weinberg <- function(fname){
    hwe <- read.table(fname,header=T)
    # hwe <- hwe[hwe$TEST=="UNAFF",]
    hist(hwe$P)
    return(hwe)
}


minor_allele_frequency <- function(fname,taille=50){
    maf <- list(data=NA,x=NA,y=NA)
    maf$data <- read.table(fname,header=T)
    maf$data$MAF[is.na(maf$data$MAF)] = 0
    maf$x <- seq(0,1, length.out=taille)
    maf$y <- c()
    for (i in 1:taille){
	maf$y <- c(sum(maf$data$MAF<=maf$x[i])/length(maf$data$MAF),maf$y)
    }
    plot(maf$y~maf$x, type='l', col='red', main="MAF")
    abline(v=0.05)
    return(maf)
}

if ( !is.null(opt$save) ) { pdf(file=ifelse(FALSE, opt$save, opt$save)) }

if ( !is.null(opt$lmiss) ) { lmiss_call_rate(opt$lmiss) }

if ( !is.null(opt$imiss) ) { imiss_call_rate(opt$imiss) }

if ( !is.null(opt$hardy) ) { hardy_weinberg(opt$hardy) }

if ( !is.null(opt$maf) ) { minor_allele_frequency(opt$maf) }

if ( !is.null(opt$save) ) { dev.off() }


