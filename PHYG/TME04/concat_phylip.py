#!env python

import argparse
import warnings
import operator

import numpy as np

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-if', '--input-file-names', dest='ifname', help='input file name list',required=True, nargs='+')
    parser.add_argument('-of', '--output-file-name', dest='ofname', help='output file name',required=True)
    return parser.parse_args()


def load_phylip_to_dict(fname,seqRec_dict={}):
    align = AlignIO.read(fname, "phylip")
    for seq_record in list(align) :
        key = seq_record.name
        if seqRec_dict.has_key(key) :
            raise ValueError("Conflict Sequences Record Key : [ %s - %s]"%(fname,seq_record.name),UserWarning)
        seqRec_dict[key] = seq_record
    return seqRec_dict


def write_dict_to_phylip(ofname,seqRec_dict):
    seqRec_list = [ rec for key,rec in seqRec_dict.iteritems() ]
    with open(ofname, 'w') as handler :
        AlignIO.write(MultipleSeqAlignment(seqRec_list), handler, "phylip")

        
def concat_to_dict(phylip_dict_list, concat_dict={}) :

    species_set = set()
    for phylip_dict in phylip_dict_list :
        for specie in phylip_dict.keys() :
            species_set.add(specie)

    concat_count = {}
    for specie in species_set :
        concat_dict[specie] = SeqRecord(Seq("", IUPAC.protein),
                                        id=specie, name=specie,
                                        description=specie)
        concat_count[specie] = 0

    def align_len_check(seqRec_dict):
        pstr = ""
        ok = 0
        rflag = False
        for k,v in seqRec_dict.iteritems():
            if ok < len(v.seq) :
                ok = len(v.seq) 
            elif ok > len(v.seq) :
                rflag = True
            pstr += "%s %d "%(k,len(v.seq))
        if rflag :
            raise ValueError(pstr)
        
    for seqRec_dict in phylip_dict_list :
        for ks,v in concat_dict.iteritems() :
            for k in ks.split("|") :
                if seqRec_dict.has_key(k) :
                    cid = concat_dict[ks].id
                    concat_dict[ks] += seqRec_dict[k]
                    concat_dict[ks].id = cid
                    concat_dict[ks].name = cid
                    concat_dict[ks].description = cid
                    concat_count[ks] += 1
                    break
                
        cmax = max(concat_count.iteritems(), key=operator.itemgetter(1))[0]
        # print concat_count
        for k,c in concat_count.iteritems():
            if concat_count[cmax] > c :
                cid = concat_dict[k].id
                concat_dict[k] += SeqRecord(Seq("-"*len(seqRec_dict[cmax].seq), IUPAC.protein),id=cid, name=cid, description=cid)
                concat_dict[k].id = cid
                concat_dict[k].name = cid
                concat_dict[k].description = cid
                concat_count[k] += 1

        align_len_check(concat_dict)
            
    return concat_dict

                
if __name__ == "__main__":
    args = interface_standard()
    phylip_dict_list = [ load_phylip_to_dict(ifname, {}) for ifname in args.ifname ]
    concat_dict = concat_to_dict(phylip_dict_list)
    write_dict_to_phylip(args.ofname, concat_dict)    
        
