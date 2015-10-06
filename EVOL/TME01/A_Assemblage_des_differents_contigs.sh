#!env bash

# include ncbi utils in PATH
export PATH=/usr/local/ncbi-blast-2.2.29+/bin/:$PATH

#### A Q.1 
# establish DB
# makeblastdb -in contigs.fasta -input_type 'fasta' -out 'contigsDB' -dbtype 'nucl'

# sequences alignment to DB
# blastn -db 'contigsDB' -query contigs.fasta -strand 'both' -out contigs.blastn -outfmt '7'

# expose assemblage possible
python assemblage_contigs_blastn_table.py contigs.blastn

# 
# blastp
