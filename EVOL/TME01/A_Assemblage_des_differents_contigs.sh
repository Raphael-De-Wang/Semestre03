#!env bash

# include ncbi utils in PATH
export PATH=/usr/local/ncbi-blast-2.2.29+/bin/:$PATH

#### A Q.1 
# establish DB
# makeblastdb -in contigs.fasta -input_type 'fasta' -out 'contigsDB' -dbtype 'nucl'

# sequences alignment to DB
# blastn -db 'contigsDB' -query contigs.fasta -strand 'both' -out contigs.blastn -outfmt '7'
## blastn -db 'contigsDB' -query contigs.fasta -strand 'both' -out contigs_alignment.blastn 

# expose assemblage possible
# python assemblage_contigs_blastn_table.py contigs.blastn

# 
python extract_sequence_start_end.py contigs.fasta 40 contigs_head_tail.fasta
makeblastdb -in contigs_head_tail.fasta -input_type 'fasta' -out 'contigsHeadTailDB' -dbtype 'nucl'
blastn -db 'contigsHeadTailDB' -query contigs_head_tail.fasta -strand 'both' -out contigs_head_tail.blastn -outfmt '6'

# combine contigs
# python combine_contigs.py -in contigs.fasta -on combine_contigs.fasta -o 8 3 2 7 4 5 -s 0 0 0 0 1 0 -L 143 98 127 100 158
# python combine_contigs.py -in contigs.fasta -on combine_contigs.fasta -o 6 1 2 7 4 5 -s 1 1 0 0 1 0 -L 122 98 127 100 158

# blastp
