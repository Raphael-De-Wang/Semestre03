#!env python

import argparse
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

# gb_list = [ "seq_EColi.gb", "seq_Salm.gb", "seq_Staphy.gb" ]
# fasta_list = [ "seq_EColi.fasta", "seq_Salm.fasta", "seq_Staphy.fasta" ]

def interface_standard():
    parser = argparse.ArgumentParser("Simple command for convert geneBank to aa fasta.")
    parser.add_argument('-gb', '--geneBank-file-name', dest='ingb', help='input geneBank file name',required=True)
    parser.add_argument('-f','--fasta-file-name', dest='infasta',help='input fasta file name',required=True)
    parser.add_argument('-o','--out-fasta-file-name', dest='outfasta',help='output fasta file name',required=True)
    return parser.parse_args()
    

def load_fasta_to_list(fname) :
    seq_list = []
    for seq_record in SeqIO.parse(fname, "fasta") :
        seq_list.append(seq_record)
    return seq_list


def load_features_to_list(fname, seq_list) :
    def nucl2aa(seq):
        return Seq.translate(seq)
    features = []
    for ind, gb_record in enumerate(SeqIO.parse(open(fname), "genbank")) :
        for f in gb_record.features:
            if cmp(f.type,"CDS") == 0 :
                seqRec = f.extract(seq_list[ind])
                seq = seqRec.seq.translate()
                try :
                    db_xref_dict = { record.split(":")[0] : record.split(":")[1] for record in f.qualifiers["db_xref"] }
                except KeyError :
                    db_xref_dict = {}
                desc= "".join([ "|%s:%s"%(key,value) for key, value in db_xref_dict.iteritems() ])
                sid = ""
                if db_xref_dict.has_key("GeneID") :
                    sid = "GeneID|" + db_xref_dict["GeneID"]
                elif db_xref_dict.has_key("GI") :
                    sid = ( "GI|" + db_xref_dict["GI"] )
                elif db_xref_dict.has_key("protein_id") :
                    sid = ( "protein_id|" + db_xref_dict["protein_id"] )
                elif db_xref_dict.has_key("PSEUDO") :
                    sid = ( "PSEUDO|" + db_xref_dict["PSEUDO"] )
                elif f.qualifiers.has_key("locus_tag") :
                    sid = ( "locus_tag|" + "".join(f.qualifiers["locus_tag"]) )
                if len(sid) == 0 :
                    warnings.warn("Feature without identity. [%s]"(fname))
                    print f
                    continue
                snm = sid
                sr = SeqRecord(seq, id=sid, name=snm, description=desc)
                if f.qualifiers.has_key('db_xref') :
                    sr.dbxrefsc = f.qualifiers['db_xref']
                features.append(sr)
    return features


def write_seq_list_to_fasta(fname, seq_list):
    SeqIO.write(seq_list, fname, "fasta")


if __name__ == "__main__":
    args = interface_standard()
    seq_list = load_fasta_to_list(args.infasta)
    features = load_features_to_list(args.ingb,seq_list)
    write_seq_list_to_fasta(args.outfasta, features)
