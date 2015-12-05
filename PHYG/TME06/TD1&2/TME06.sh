#!env bash

# currpath=`pwd`
# cd programs/
# tar -xf mga_x86_64.tar.gz 
# cd $currpath

# Ex 02

mga="programs/mga_linux_ia64"

for genom in `ls genomes/*fna` ; do 
    echo "$mga $genom > ${genom/fna/mga}"
done


# Ex 03

for genom in `ls genomes/*fna` ; do 
    python mgaParser.py -g $genom -i ${genom/fna/mga} -o ${genom/fna/fasta}
done

# Ex 04

python ex04.py -if `ls genomes/*fasta` -p ex04.pdf

# Ex 05

python ex05.py -if `ls genomes/*fasta` -p ex05.pdf

# Ex 06
echo "pretending to read guide manual..." 

# Ex 07
# currpath=`pwd`
# cd programs/cd-hit-auxtools/
# make
# cd $currpath

# Ex 08
cat genomes/*fasta > all.fasta

# Ex 09
cd_hit="programs/cdhit-master/cd-hit"
echo "utlisation of cd-hit refer to page 4 to 6 in cd-hit_user_s_guide.pdf"
$cd_hit -i all.fasta -o analysis.out -c 0.7
# python ex09.py -f all.fasta


