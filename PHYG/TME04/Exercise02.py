#env python

import argparse
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-if', '--input-file-name', dest='ifname', help='input fasta name',required=True, nargs='+')
    parser.add_argument('-is', '--input-species-list', dest='isname', help='input species list file name',required=True)
    parser.add_argument('-of', '--output-file-name', dest='ofname', help='output fasta file name',required=True)
    return parser.parse_args()


# Q2.1
def load_fasta_to_dict(fname,seqRec_dict):
    for seq_record in SeqIO.parse(fname, "fasta"):
        key = seq_record.name + '_' + fname.split(".")[0]
        if seqRec_dict.has_key(key) :
            warnings.warn("Conflict Sequences Record Key : [ %s - %s]"%(fname,seq_record.name),UserWarning)
        seqRec_dict[key] = seq_record
    return seqRec_dict


# Q2.2
def load_species_list(fname):
    species_set = set()
    with open(fname) as f:
        for line in f :
            for specie in line.strip().split(";") :
                if len(specie) > 0 : 
                    species_set.add(specie)
    return list(species_set)
        

def extract_species_seqRec_list(seqRec_dict, species_list):
    seqRec_list = []
    for sid,rec in seqRec_dict.iteritems():
        for ssid in sid.split("_") :
            if ssid in species_list :
                rec = seqRec_dict[sid]
                id = "|".join(rec.id.split("_")[::-1])
                rec.id = id
                rec.name = id
                rec.description = id
                seqRec_list.append(rec)
    return seqRec_list 

                
def write_seqRec_list_to_fasta(seqRec_list,fname):
    SeqIO.write(seqRec_list, fname, "fasta")
                  
                  
if __name__ == "__main__":
    args = interface_standard()
    seqRec_dict = {}
    for ifname in args.ifname :
        load_fasta_to_dict(ifname,seqRec_dict)
    species_list = load_species_list(args.isname)
    seqRec_list = extract_species_seqRec_list(seqRec_dict, species_list)
    write_seqRec_list_to_fasta(seqRec_list,args.ofname)