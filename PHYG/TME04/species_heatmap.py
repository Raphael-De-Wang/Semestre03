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
    parser.add_argument('-if', '--input-file-name', dest='ifname', help='input file name list',required=True, nargs='+')
    parser.add_argument('-is', '--input-species-list', dest='isname', help='input species file name',required=True)
    parser.add_argument('-of', '--heatmap-file-name', dest='ofname', help='heatmap plot file name',required=True)
    parser.add_argument('--no-heatmap', dest='heatmap', help='active heatmap plot, and return plot address', action='store_false')
    parser.add_argument('--filter', dest='filter', help='active filter mode, weak features will be elimated by given thresholds of data lost percentage in pfams and species', nargs='+')
    parser.add_argument('--concat', dest='concat', help='active concat mode, filtered alignments will be concat to named file')
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

def load_species_list(fname):
    species_label = []
    with open(fname) as f:
        for line in f :
            species = line.strip().split(";")
            if '' in species :
                species.remove('')
            species_label.append(species)
    return species_label


def build_heatmap_mtx(keys_list, species_list) :
    heatmap = []
    for keys in keys_list :
        hv = []
        for specie in species_list :
            scount = 0
            for s in specie :
                if s in keys :
                    scount += 1
                    break
            hv.append(scount)
        heatmap.append(hv)
    return heatmap


def plot_heatmap(x,y,z,ofname=None,eli_x=[],eli_y=[]):
    import plotly.plotly as py
    import plotly.graph_objs as go
    py.sign_in("WangJinxin", "qy753qduyy")
    if not isinstance(eli_x,list) or not isinstance(eli_y,list) :
        raise ValueError("Incorrect Data Type")
    x_ind = [ i for i in range(len(x)) if i not in eli_x ]
    y_ind = [ i for i in range(len(y)) if i not in eli_y ]
    # data = [ go.Heatmap(x=x,y=y,z=z) ]
    data = [ go.Heatmap(x=np.array(x)[x_ind],y=np.array(y)[y_ind],z=np.array(z)[y_ind,:][:,x_ind]) ]
    plot_url = py.plot(data, filename=ofname, auto_open=False)
    print plot_url+".embed"


def remove_weak_axis(heatmap, axis=0, seuil=0.5) :
    weak_index = np.where(np.sum(heatmap,axis=axis)*1./np.shape(heatmap)[axis] < seuil)
    return weak_index[0].tolist()


def concat_to_dict(pfams, species_list, eli_pfams=[], eli_species=[], concat_dict={}) :
    if not isinstance(eli_pfams,list) or not isinstance(eli_species,list) :
        raise ValueError("Incorrect Data Type")
    
    pfams_ind = [ i for i in range(len(pfams)) if i not in eli_pfams ]
    species_ind = [ i for i in range(len(species)) if i not in eli_species ]

    concat_count = {}
    for ind in species_ind :
        cid = species_list[ind][0]
        concat_dict["|".join(species_list[ind])] = SeqRecord(Seq("", IUPAC.protein),
                                                             id=cid, name=cid,
                                                             description=cid)
        concat_count["|".join(species_list[ind])] = 0

    def simplify_seqRec_dict(seqRec_dict,seqRec_sim_dict={}):
        for k,v in seqRec_dict.iteritems():
            seqRec_sim_dict[k.split("|")[0]] = v
        return seqRec_sim_dict

    def align_len_check(pfam, seqRec_dict):
        pstr = pfam + " : "
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
        
    for pfam in np.array(pfams)[pfams_ind] :
        fname = pfam + ".phylip"
        seqRec_dict = {}
        seqRec_sim_dict = {}
        load_phylip_to_dict(fname,seqRec_dict)
        simplify_seqRec_dict(seqRec_dict,seqRec_sim_dict)
        for ks,v in concat_dict.iteritems() :
            for k in ks.split("|") :
                if seqRec_sim_dict.has_key(k) :
                    cid = concat_dict[ks].id
                    concat_dict[ks] += seqRec_sim_dict[k]
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
                concat_dict[k] += SeqRecord(Seq("-"*len(seqRec_sim_dict[cmax].seq), IUPAC.protein),id=cid, name=cid, description=cid)
                concat_dict[k].id = cid
                concat_dict[k].name = cid
                concat_dict[k].description = cid
                concat_count[k] += 1

        align_len_check(pfam, concat_dict)
            
    return concat_dict

                
if __name__ == "__main__":
    args = interface_standard()
    species_list = load_species_list(args.isname)
    keys_list = []
    eli_pfams = []
    eli_species = []
    for ifname in args.ifname :
        seqRec_dict = {}
        load_phylip_to_dict(ifname,seqRec_dict)
        keys = [ key.split("|")[0] for key in seqRec_dict.keys() ]
        keys_list.append(keys)
    heatmap = build_heatmap_mtx(keys_list, species_list)
    pfams   = [ fname.split('.')[0] for fname in args.ifname ]
    species = [ species[0] for species in species_list ]
    heatmap = np.array(heatmap)
    
    if args.filter :
        pfams_seuil =   float(args.filter[0])
        species_seuil = float(args.filter[1])
        eli_pfams   = remove_weak_axis(heatmap.T, axis=0, seuil=pfams_seuil)
        eli_species = remove_weak_axis(heatmap.T, axis=1, seuil=species_seuil)
        
    if args.heatmap :
        plot_heatmap(pfams, species, heatmap.T, ofname=args.ofname, eli_x=eli_pfams, eli_y=eli_species)

    if args.concat :
        concat_dict = concat_to_dict(pfams, species_list, eli_pfams=eli_pfams, eli_species=eli_species)
        write_dict_to_phylip(args.concat, concat_dict)    
        
