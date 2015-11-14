#!env bash

# include ncbi utils in PATH
export PATH=/usr/local/ncbi-blast-2.2.29+/bin/:$PATH

# Escherichia_coli.fasta  
# Salmonella_enterica.fasta  
# Staphylococcus_aureus.fasta

#### Q.A : Recuperer les génomes à étudier
# python gb2aa.py -gb seq_EColi.gb -f seq_EColi.fasta -o seq_EColi_aa.fasta 
# python gb2aa.py -gb seq_Salm.gb -f seq_Salm.fasta -o seq_Salm_aa.fasta 
# python gb2aa.py -gb seq_Staphy.gb -f seq_Staphy.fasta -o seq_Staphy_aa.fasta 


#### Q.B : Construction des familles de genes

# makeblastdb -in Escherichia_coli_prot.fasta -input_type 'fasta' -out 'EscherichiaColiDB' -dbtype 'prot'
# makeblastdb -in Salmonella_enterica_prot.fasta -input_type 'fasta' -out 'SalmonellaEntericaDB' -dbtype 'prot'
# makeblastdb -in Staphylococcus_aureus_prot.fasta -input_type 'fasta' -out 'StaphylococcusAureusDB' -dbtype 'prot'

# blastp -db 'EscherichiaColiDB' -query Salmonella_enterica_prot.fasta  -out Sa_vs_Es.blastn -outfmt '6'
# blastp -db 'EscherichiaColiDB' -query Staphylococcus_aureus_prot.fasta  -out St_vs_Es.blastn -outfmt '6'

# blastp -db 'SalmonellaEntericaDB' -query Escherichia_coli_prot.fasta  -out Es_vs_Sa.blastn -outfmt '6'
# blastp -db 'SalmonellaEntericaDB' -query Staphylococcus_aureus_prot.fasta  -out St_vs_Sa.blastn -outfmt '6'

# blastp -db 'StaphylococcusAureusDB' -query Escherichia_coli_prot.fasta  -out Es_vs_St.blastn -outfmt '6'
# blastp -db 'StaphylococcusAureusDB' -query Salmonella_enterica_prot.fasta  -out Sa_vs_St.blastn -outfmt '6'

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

# silix -i $iden -r $overlap Escherichia_coli_prot.fasta Sa_vs_Es_best_hits.blastn > EColi_Sa_vs_Es_best_hits.silix
# silix -i $iden -r $overlap Escherichia_coli_prot.fasta Es_vs_Sa_best_hits.blastn > EColi_Es_vs_Sa_best_hits.silix

# silix -i $iden -r $overlap Escherichia_coli_prot.fasta St_vs_Es_best_hits.blastn > EColi_St_vs_Es_best_hits.silix
# silix -i $iden -r $overlap Escherichia_coli_prot.fasta Es_vs_St_best_hits.blastn > EColi_Es_vs_St_best_hits.silix

# silix -i $iden -r $overlap Staphylococcus_aureus_prot.fasta Sa_vs_Es_best_hits.blastn> Staphy_Sa_vs_Es_best_hits.silix
# silix -i $iden -r $overlap Staphylococcus_aureus_prot.fasta Es_vs_Sa_best_hits.blastn> Staphy_Es_vs_Sa_best_hits.silix

# silix -i $iden -r $overlap Staphylococcus_aureus_prot.fasta Sa_vs_St_best_hits.blastn> Staphy_Sa_vs_St_best_hits.silix
# silix -i $iden -r $overlap Staphylococcus_aureus_prot.fasta St_vs_Sa_best_hits.blastn> Staphy_St_vs_Sa_best_hits.silix

# silix -i $iden -r $overlap Salmonella_enterica_prot.fasta St_vs_Sa_best_hits.blastn  > Salm_St_vs_Sa_best_hits.silix
# silix -i $iden -r $overlap Salmonella_enterica_prot.fasta Sa_vs_St_best_hits.blastn  > Salm_Sa_vs_St_best_hits.silix

# silix -i $iden -r $overlap Salmonella_enterica_prot.fasta St_vs_Es_best_hits.blastn  > Salm_St_vs_Es_best_hits.silix
# silix -i $iden -r $overlap Salmonella_enterica_prot.fasta Es_vs_St_best_hits.blastn  > Salm_Es_vs_St_best_hits.silix

wc -l *best_hits.blastn
wc -l *silix

cat *best_hits.blastn > all_best_hits.blastn
cat *prot.fasta > all_prot.fasta

silix all_prot.fasta all_best_hits.blastn > all_fam.silix

#### Q.C : Analyse

python dist_hist.py 

python venn_hist.py -f all_fam.silix  -o venn3.png
