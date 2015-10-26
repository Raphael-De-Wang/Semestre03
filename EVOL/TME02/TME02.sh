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

makeblastdb -in Escherichia_coli_prot.fasta -input_type 'fasta' -out 'EscherichiaColiDB' -dbtype 'prot'
makeblastdb -in Salmonella_enterica_prot.fasta -input_type 'fasta' -out 'SalmonellaEntericaDB' -dbtype 'prot'
makeblastdb -in Staphylococcus_aureus_prot.fasta -input_type 'fasta' -out 'StaphylococcusAureusDB' -dbtype 'prot'

blastp -db 'EscherichiaColiDB' -query Salmonella_enterica_prot.fasta  -out Sa_vs_Es.blastn -outfmt '6'
blastp -db 'EscherichiaColiDB' -query Staphylococcus_aureus_prot.fasta  -out St_vs_Es.blastn -outfmt '6'

blastp -db 'SalmonellaEntericaDB' -query Escherichia_coli_prot.fasta  -out Es_vs_Sa.blastn -outfmt '6'
blastp -db 'SalmonellaEntericaDB' -query Staphylococcus_aureus_prot.fasta  -out St_vs_Sa.blastn -outfmt '6'

blastp -db 'StaphylococcusAureusDB' -query Escherichia_coli_prot.fasta  -out Es_vs_St.blastn -outfmt '6'
blastp -db 'StaphylococcusAureusDB' -query Salmonella_enterica_prot.fasta  -out Sa_vs_St.blastn -outfmt '6'

#### python select_best_blast_hits.py ####

# silix -i Escherichia_coli.fasta Sa_vs_Es.blastn
# silix -i Escherichia_coli.fasta Es_vs_Sa.blastn

# silix -i Escherichia_coli.fasta St_vs_Es.blastn
# silix -i Escherichia_coli.fasta Es_vs_St.blastn

# silix -i Salmonella_enterica.fasta Sa_vs_Es.blastn
# silix -i Salmonella_enterica.fasta Es_vs_Sa.blastn

# silix -i Salmonella_enterica.fasta Sa_vs_St.blastn
# silix -i Salmonella_enterica.fasta St_vs_Sa.blastn

# silix -i Staphylococcus_aureus.fasta St_vs_Sa.blastn
# silix -i Staphylococcus_aureus.fasta Sa_vs_St.blastn

# silix -i Staphylococcus_aureus.fasta St_vs_Es.blastn
# silix -i Staphylococcus_aureus.fasta Es_vs_St.blastn

#### Q.C : Analyse

# python dist_hist.py

# python venn_hist.py
