#!env bash

# include ncbi utils in PATH
export PATH=/usr/local/ncbi-blast-2.2.29+/bin/:$PATH

# Escherichia_coli.fasta  
# Salmonella_enterica.fasta  
# Staphylococcus_aureus.fasta

#### Q.A : Recuperer les génomes à étudier
python gb2aa.py -gb seq_EColi.gb -f seq_EColi.fasta -o seq_EColi_aa.fasta 
python gb2aa.py -gb seq_Salm.gb -f seq_Salm.fasta -o seq_Salm_aa.fasta 
python gb2aa.py -gb seq_Staphy.gb -f seq_Staphy.fasta -o seq_Staphy_aa.fasta 


#### Q.B : Construction des familles de genes

makeblastdb -in seq_EColi_aa.fasta -input_type 'fasta' -out 'EscherichiaColiDB' -dbtype 'prot'
makeblastdb -in seq_Salm_aa.fasta -input_type 'fasta' -out 'SalmonellaEntericaDB' -dbtype 'prot'
makeblastdb -in seq_Staphy_aa.fasta -input_type 'fasta' -out 'StaphylococcusAureusDB' -dbtype 'prot'

blastp -db 'EscherichiaColiDB' -query seq_Salm_aa.fasta  -out Sa_vs_Es.blastn -outfmt '6'
blastp -db 'EscherichiaColiDB' -query seq_Staphy_aa.fasta  -out St_vs_Es.blastn -outfmt '6'

blastp -db 'SalmonellaEntericaDB' -query seq_EColi_aa.fasta  -out Es_vs_Sa.blastn -outfmt '6'
blastp -db 'SalmonellaEntericaDB' -query seq_Staphy_aa.fasta  -out St_vs_Sa.blastn -outfmt '6'

blastp -db 'StaphylococcusAureusDB' -query seq_EColi_aa.fasta  -out Es_vs_St.blastn -outfmt '6'
blastp -db 'StaphylococcusAureusDB' -query seq_Salm_aa.fasta  -out Sa_vs_St.blastn -outfmt '6'

#### Dresser la liste des Reciprocal Best Hits entre les 3 paires de génomes.

iden=0
evalue=1.0

python select_best_blast_hits.py -f Sa_vs_Es.blastn -o Sa_vs_Es_best_hits.blastn -i $iden -e $evalue
python select_best_blast_hits.py -f St_vs_Es.blastn -o St_vs_Es_best_hits.blastn -i $iden -e $evalue

python select_best_blast_hits.py -f Es_vs_Sa.blastn -o Es_vs_Sa_best_hits.blastn -i $iden -e $evalue 
python select_best_blast_hits.py -f St_vs_Sa.blastn -o St_vs_Sa_best_hits.blastn -i $iden -e $evalue 

python select_best_blast_hits.py -f Es_vs_St.blastn -o Es_vs_St_best_hits.blastn -i $iden -e $evalue 
python select_best_blast_hits.py -f Sa_vs_St.blastn -o Sa_vs_St_best_hits.blastn -i $iden -e $evalue 


#### Dresser la liste des familles de gènes homologues partagés par ces 3 génomes.

iden=0.35
overlap=0.80

silix -i $iden -r $overlap seq_EColi_aa.fasta Sa_vs_Es_best_hits.blastn > EColi_Sa_vs_Es_best_hits.silix
silix -i $iden -r $overlap seq_EColi_aa.fasta Es_vs_Sa_best_hits.blastn > EColi_Es_vs_Sa_best_hits.silix

silix -i $iden -r $overlap seq_EColi_aa.fasta St_vs_Es_best_hits.blastn > EColi_St_vs_Es_best_hits.silix
silix -i $iden -r $overlap seq_EColi_aa.fasta Es_vs_St_best_hits.blastn > EColi_Es_vs_St_best_hits.silix

silix -i $iden -r $overlap seq_Staphy_aa.fasta Sa_vs_Es_best_hits.blastn> Staphy_Sa_vs_Es_best_hits.silix
silix -i $iden -r $overlap seq_Staphy_aa.fasta Es_vs_Sa_best_hits.blastn> Staphy_Es_vs_Sa_best_hits.silix

silix -i $iden -r $overlap seq_Staphy_aa.fasta Sa_vs_St_best_hits.blastn> Staphy_Sa_vs_St_best_hits.silix
silix -i $iden -r $overlap seq_Staphy_aa.fasta St_vs_Sa_best_hits.blastn> Staphy_St_vs_Sa_best_hits.silix

silix -i $iden -r $overlap seq_Salm_aa.fasta St_vs_Sa_best_hits.blastn  > Salm_St_vs_Sa_best_hits.silix
silix -i $iden -r $overlap seq_Salm_aa.fasta Sa_vs_St_best_hits.blastn  > Salm_Sa_vs_St_best_hits.silix

silix -i $iden -r $overlap seq_Salm_aa.fasta St_vs_Es_best_hits.blastn  > Salm_St_vs_Es_best_hits.silix
silix -i $iden -r $overlap seq_Salm_aa.fasta Es_vs_St_best_hits.blastn  > Salm_Es_vs_St_best_hits.silix

# wc -l *silix

#### Q.C : Analyse

python dist_hist.py

python venn_hist.py
