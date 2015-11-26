#!env python

import argparse
import warnings

from Bio import AlignIO

import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-if', '--input-file-name', dest='ifname', help='input file name list',required=True, nargs='+')
    parser.add_argument('-is', '--input-species-list', dest='isname', help='input species file name',required=True)
    parser.add_argument('-of', '--heatmap-file-name', dest='ofname', help='heatmap plot file name',required=True)
    parser.add_argument('--filter', dest='filter', help='active filter mode, weak features will be elimated', action='store_true')
    return parser.parse_args()


def load_phylip_to_dict(fname,seqRec_dict):
    align = AlignIO.read(fname, "phylip")
    for seq_record in list(align) :
        key = seq_record.name
        if seqRec_dict.has_key(key) :
            warnings.warn("Conflict Sequences Record Key : [ %s - %s]"%(fname,seq_record.name),UserWarning)
        seqRec_dict[key] = seq_record
    return seqRec_dict


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
            hv.append(scount)
        heatmap.append(hv)
    return heatmap


def plot_heatmap(x,y,z,ofname=None,eli_x=[],eli_y=[]):
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
    # plot_heatmap([ species[0] for species in species_list ], [ fname.split('.')[0] for fname in args.ifname ],heatmap,ofname=args.ofname)
    heatmap = np.array(heatmap)
    if args.filter : 
        eli_pfams   = remove_weak_axis(heatmap.T, axis=0, seuil=0.9)
        eli_species = remove_weak_axis(heatmap.T, axis=1, seuil=0.8)
    plot_heatmap(pfams, species, heatmap.T, ofname=args.ofname, eli_x=eli_pfams, eli_y=eli_species)

