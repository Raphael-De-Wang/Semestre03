library(GenABEL)


#==============================================================================
# Prepare Data                                                                          
#==============================================================================

workpath <- "/users/Etu9/3404759/Workspace/Semestre03/GENOM/"

concatPhenoCovar <- function () {
    pheno_1 = read.table(paste(workpath,"Data/phenotype_plink.txt",sep=''),header=T)
    pheno_2 = read.table(paste(workpath,"Data/covar_plink.txt",sep=''),header=T)
    pheno_file = matrix(ncol=6,nrow=3175)
    pheno_file[,1] = pheno_1[,1]
    pheno_file[,2] = as.numeric(pheno_1[,3])
    pheno_file[,3] = as.numeric(pheno_2[,4])
    pheno_file[,4] = as.numeric(pheno_2[,3])-1
    pheno_file[,5] = as.numeric(pheno_2[,4])
    pheno_file[,5] = as.numeric(pheno_2[,4])
    pheno_file[,6] = as.numeric(pheno_2[,5])
    pheno_file[,1] = as.character(pheno_1[,1])
    colnames(pheno_file)=c("id","ApoA1","HDL_CHOL","sex","age","BMI")
    write.table(pheno_file,paste(workpath,"Data/pheno_file.txt",sep=''),row.names = F)
    return(pheno_file)
}
 
relatedness <- function () {
    d = read.table("../Data/QC/chr13_relatedness.genome", header=T)
    par(pch=16)
    with(d,plot(d$Z0,d$Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
    with(subset(d,d$RT=="FS") , points(d$Z0,d$Z1,col=3))
    with(subset(d,d$RT=="HS") , points(d$Z0,d$Z1,col="darkorange"))
    with(subset(d,d$RT=="OT") , points(d$Z0,d$Z1,col=4))
    with(subset(d,d$RT=="PO") , points(d$Z0,d$Z1,col=2))
    with(subset(d,d$RT=="UN") , points(d$Z0,d$Z1,col=1))
    with(subset(d,d$RT=="OT") , points(d$Z0,d$Z1,col=4))
    with(subset(d,d$RT=="HS") , points(d$Z0,d$Z1,col="darkorange"))

    legend(1,1, xjust=1, yjust=1, legend=levels(d$RT), pch=16, col=c(3,"darkorange",4,2,1))
    library(ggplot2)
    qplot(d$Z0,d$Z1,data=d,colour=d$RT)
    scale_colour_manual(value=c(3,"darkorange",4,2,1))
}

pedfile <- "Data/QC/chr13_without_missing.txt"
# data_test = convert.snp.ped(paste(workpath,pedfile,sep=''),paste(workpath,"Data/chr13.map",sep=''),paste(workpath,"Data/QC/chr13_geno_file.txt",sep=''))
concatPhenoCovar()
gwaa.data <- load.gwaa.data(phenofile = paste(workpath,"Data/pheno_file.txt",sep=''),genofile = paste(workpath,"Data/QC/chr13_geno_file.txt",sep=''))
check.marker(gwaa.data)
