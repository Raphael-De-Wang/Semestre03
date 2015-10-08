#!env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def load_input_fasta_to_list(fname):
    seqRec_list = []
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqRec_list.append(seq_record)
    return seqRec_list

def save_combined_contigs_to_fasta(fname,record_list):
    SeqIO.write(record_list, fname, "fasta")
    
ifname = sys.argv[1]
seuil = int(sys.argv[2])
ofname = sys.argv[3]
seqRec_list = load_input_fasta_to_list(ifname)
record_list = []

for seq_record in seqRec_list :
    record = seq_record[:seuil]
    record.id += " head %d"%seuil
    record.name += " head %d"%seuil
    record_list.append(record)

    record_comp = record.reverse_complement()
    record_comp.id = record.id + "compl"
    record_comp.name = record.name + "compl"
    record_list.append(record_comp)

    record = seq_record[-seuil:]
    record.id += " tail %d"%seuil
    record.name += " tail %d"%seuil
    record_list.append(record)

    record_comp = record.reverse_complement()
    record_comp.id = record.id + "compl"
    record_comp.name = record.name + "compl"
    record_list.append(record_comp)
    
save_combined_contigs_to_fasta(ofname,record_list)