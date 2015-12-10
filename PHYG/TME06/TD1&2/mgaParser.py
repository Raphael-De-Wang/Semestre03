#!env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein


class mga_record (object) :
    def __init__(self, record_line,genome_ind) :
        self.record = record_line.split()
        try :
            self.genome_ind= genome_ind
            self.gene_id   = self.record[0]
            self.start_pos = int(self.record[1])
            self.end_pos   = int(self.record[2])
            self.strand    = self.record[3]
            self.frame     = int(self.record[4])
            self.comp_part = int(self.record[5])
            self.gene_score= float(self.record[6])
            self.used_model= self.record[7]
        except :
            print record_line
            raise
        
        if self.record[8] == "-" :
            self.rbs_start = self.record[8]
        else :
            self.rbs_start = int(self.record[8])
            
        if self.record[9] == "-" :            
            self.rbs_end   = self.record[9]
        else :
            self.rbs_start = int(self.record[9])
            
        if self.record[10] == "-" :
            self.rbs_score = self.record[10]
        else :
            self.rbs_score = float(self.record[10])
            
    def __str__(self):
        return " ".join(self.record)

def interface_standard():
    parser = argparse.ArgumentParser("MetaGeneAnnotator Parser Script")
    parser.add_argument('-i', '--input',  dest='inputFile', help='Gene record list',required=True)
    parser.add_argument('-o', '--output', dest='outputFile', help='fasta output',required=True)
    parser.add_argument('-g', '--genome', dest='genomeFile', help='MetaGeneAnnotator output',required=True)
    return parser.parse_args()


def mga_parser(genomeFile,inputFile) :

    gene_list = []
    gene_records = []
    genome_list = [ seq_record for seq_record in SeqIO.parse(genomeFile, "fasta") ]

    gfname = genomeFile.split("/")[-1]
    genom_name = gfname.split("_")[-2]
    if len(genom_name) < 3 :
        genom_name = gfname.split("_")[-3]
        
    with open(inputFile) as handler :
        comment_count = 0
        for record_line in handler :
            if cmp(record_line[0],'#') == 0 :
                comment_count += 1
                continue
            gene_records.append(mga_record(record_line,(comment_count/3)-1))
    
    for mga in gene_records :
        seq = genome_list[mga.genome_ind].seq
        if mga.strand == "+" :
            #seq = seq.reverse_complement()
            aa = seq[mga.start_pos - 1 : mga.end_pos].translate()
        else:
            aa = seq[mga.start_pos - 1 : mga.end_pos].reverse_complement().translate()
        seq_name = mga.gene_id + "_" + genom_name
        seq_rec = SeqRecord(aa, id=seq_name, name=seq_name, description="")
        gene_list.append(seq_rec)
        
    return gene_list
    

if __name__ == "__main__":
    args = interface_standard()
    gene_list = mga_parser(args.genomeFile,args.inputFile)
    SeqIO.write(gene_list, open(args.outputFile,"w"), "fasta")
	
