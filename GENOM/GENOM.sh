#!env bash

plink="/users/Etu9/3404759/Workspace/toolkits/plink_linux_1.90/plink"
prettify="/users/Etu9/3404759/Workspace/toolkits/plink_linux_1.90/prettify"

QC_analyse () {

    chr=$1

    echo "  -- QC ANALYSE -- "
    
    # $plink --ped Data/QC/${chr}_without_missing.ped --map Data/$chr.map --silent --covar Data/covar_plink.txt --pheno Data/phenotype_plink.txt  --out Data/QC/$chr --missing --freq --hardy

    # Missing genotypes, Allele frequency - MAF - Minor allele frequency, Hardy-Weinberg Equilibrium
    echo "Rscript step01/qc_analysis.R --lmiss Data/QC/$chr.lmiss --imiss Data/QC/$chr.imiss --hardy Data/QC/$chr.hwe --maf Data/QC/$chr.frq --save plots/${chr}_qc_analysis.pdf > /dev/null"

}


QC_pipeline () {
    echo "  -- QC PIPELINE -- "

}

step1() {

    for chr in "chr2" "chr5" "chr13" "chr16" ; do 

	echo "---- STEP ONE QUALITY CONTROL ----"
	python step01/check_data.py -pheno Data/phenotype_plink.txt -ped Data/$chr.ped -of Data/QC/${chr}_without_missing.ped

	QC_analyse "$chr"
	# QC_pipeline
    done
}


step1
