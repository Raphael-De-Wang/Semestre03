#!/usr/bin/env python


import sys
import Bio.PDB
import commands
import numpy as np
from copy import deepcopy


class Vertex:
    def __init__(self, S, E):
        """S and E are Atoms objects from a Bio.PDB.Structure"""
        self.S = S
        self.E = E
        self.coord = S.coord
        self.x = S.coord[0]
        self.y = S.coord[1]
        self.z = S.coord[2]
        self.resname = S.parent.resname
        self.id = S.parent.id[1]

        
    def __getitem__(self, index):
        if index == 0:
            return self.S
        elif index == 1:
            return self.E
        elif index == 2:
            return self.id
        elif index == 3:
            return self.resname
        else:
            raise IndexError("index out of range")
    
    def __repr__(self):
        #return "<Vertex %s %s>" %(self.resname, self.id)
        return "%s %s" %(self.resname, self.id)
    
    @property
    def shortname(self):
        return "%s%s" % (Bio.PDB.protein_letters_3to1[self.resname], self.id)


class Weight:
    def __init__(self, vertex1, vertex2):
        self.residues = [vertex1.resname, vertex2.resname]
        # eudlidean distances between CA - CA and CA - CB
        self.S1S2 = np.linalg.norm(vertex1.S.coord - vertex2.S.coord)
        self.S1E2 = np.linalg.norm(vertex1.S.coord - vertex2.E.coord)
        self.S2E1 = np.linalg.norm(vertex1.E.coord - vertex2.S.coord)
        self.q = self.get_dihedral(vertex1.E.coord, vertex1.S.coord,
                                   vertex2.S.coord, vertex2.E.coord)
        #self.q = Bio.PDB.calc_dihedral()

    @classmethod
    def get_dihedral(cls, a,b,c,d):
        u = np.cross(a-b, b-c)
        v = np.cross(b-c, c-d)
        return np.rad2deg(np.arccos(np.dot(u/np.linalg.norm(u),
                                           v/np.linalg.norm(v))))

    def __repr__(self):
        return "%s, %s, %s, %s" % (self.S1S2, self.S1E2, self.S2E1, self.q)

    def reverse(self):
        self.residues.reverse()
        self.S1E2, self.S2E1 = self.S2E1, self.S1E2
        # self.q stays the same !


class Edge:
    def __init__(self, vertex1, vertex2):
        self.v1 = vertex1
        self.v2 = vertex2
        self.w = Weight(vertex1, vertex2)

    def __getitem__(self, index):
        if index == 0:
            return self.v1
        elif index == 1:
            return self.v2
        else:
            raise IndexError("index out of range")

    def __repr__(self):
        return "%s <--> %s [%s]" %(self.v1, self.v2, self.w)
    
    def reverse(self):
        self.v1, self.v2 = self.v2, self.v1
        self.w.reverse()
    
    @property
    def tabrepr(self):
        """ Vertex1 x1 y1 z1 Vertex2 x2 y2 z2 """
        return "\t".join([self.v1.shortname,
                          str(self.v1.x), str(self.v1.y), str(self.v1.z),
                          self.v2.shortname,
                          str(self.v2.x), str(self.v2.y), str(self.v2.z)]) + "\n"
    
        
#### Part 1 : construction of graph
# 2CPK.pdb
# 3LCK.pdb
####


def load_PDB(pdb_fname):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure(pdb_fname.replace(".pdb",""),pdb_fname)
    return structure


def select_aa_acessible(pdb_fname, chain_name=None, Wacc=1):

    dssp_fname = pdb_fname.replace(".pdb",".dssp")

    if commands.getstatusoutput("mkdssp -i %s -o %s"%(pdb_fname,dssp_fname))[0]:
        raise RuntimeError("mkdssp failed !")

    dssp_dict, key_list = Bio.PDB.make_dssp_dict(dssp_fname)

    if chain_name :
        return [k[1][1] for k in key_list if chain_name == k[0] and dssp_dict[k][2] >= Wacc ]
    else :
        return [k[1][1] for k in key_list if dssp_dict[k][2] >= Wacc ]
    


def select_vertices(structure, accessible_aa = None):

    vertices_list = []

    if accessible_aa is None:
        is_accessible = lambda residu: True
    else:
        is_accessible = lambda residu: residu.id[1] in accessible_aa

    for residu in structure[0].get_residues() :

        if Bio.PDB.is_aa(residu.resname) and is_accessible(residu):
            atom_pair_dict = {}

            for a in residu :
                if a.name in ['CA', 'CB', 'N'] :
                    atom_pair_dict[a.name] = a 
            
            if atom_pair_dict.get("CA") != None and atom_pair_dict.get("CB") != None:
                vertices_list.append(Vertex(atom_pair_dict["CA"], atom_pair_dict["CB"]))

            elif atom_pair_dict.get("CA")!= None and atom_pair_dict.get("N")!= None:
                vertices_list.append(Vertex(atom_pair_dict["CA"], atom_pair_dict["N"]))

            else :
                print atom_pair_dict
                print residu.id[1]
                raise ValueError("Missing value in vertex.")
    return vertices_list
    

def build_graph(vertices_list, Dmax=8): 
    graph = []
    N = len(vertices_list)
    for i in range(N):
        for j in range(i+1, N):
            a = vertices_list[i].coord
            b = vertices_list[j].coord
            dist = np.linalg.norm(a-b)
            if dist <= Dmax:
                graph.append(Edge(vertices_list[i], vertices_list[j]))
    return graph


#### Part 2 : comparison of graph

# from: http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
_blosum62_aa =  ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                 "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*" ]

    
_blosum62 = np.matrix(
""" 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ; 
   -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ;
   -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ;
   -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ;
    0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ;
   -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ;
   -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ;
    0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ;
   -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ;
   -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ;
   -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ;
   -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ;
   -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ;
   -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ;
   -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ;
    1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ;
    0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ;
   -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ;
   -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ;
    0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ;
   -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ;
   -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ;
    0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ;
   -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 """)


_compatible_aa = _blosum62 >= 2


def residu_identical(edge1, edge2, compatible_aa = _compatible_aa):
    """compatible_aa: matrix of similar amino-acids. ex: BLOSUM62 >= 2"""
    # get amino-acids of each vertex in edge
    e1_v1 = Bio.PDB.protein_letters_3to1[edge1.v1[3]]
    e1_v2 = Bio.PDB.protein_letters_3to1[edge1.v2[3]]
    e2_v1 = Bio.PDB.protein_letters_3to1[edge2.v1[3]]
    e2_v2 = Bio.PDB.protein_letters_3to1[edge2.v2[3]]
    
    if compatible_aa[_blosum62_aa.index(e1_v1), _blosum62_aa.index(e2_v1)] \
       and compatible_aa[_blosum62_aa.index(e1_v2), _blosum62_aa.index(e2_v2)]:
        return True
    if compatible_aa[_blosum62_aa.index(e1_v1), _blosum62_aa.index(e2_v2)] \
       and compatible_aa[_blosum62_aa.index(e1_v2), _blosum62_aa.index(e2_v1)]:
        edge2.reverse()
        return True

    return False


def weight_identical(edge1, edge2, DSS=0.6, DSE=0.75, Dq=35):
    """edge1 and edge2 must have been aligned using residu_identical"""
    deltaSS = abs(edge1.w.S1S2 - edge2.w.S1S2)
    deltaSE = np.array([abs(edge1.w.S1E2 - edge2.w.S1E2),
                        abs(edge1.w.S2E1 - edge2.w.S2E1)])
    deltaq = abs(edge1.w.q - edge2.w.q)

    if deltaSS <= DSS and all(deltaSE <= DSE) and deltaq <= Dq:
        return True

    return False


def edges_comp(graph1, graph2):
    common_edges_g1 = {}
    common_edges_g2 = {}
    for edge1 in graph1:
        for edge2 in graph2:
            if residu_identical(edge1,edge2) and weight_identical(edge1,edge2):
                if common_edges_g1.has_key(edge1) : 
                    common_edges_g1[edge1].append(edge2)
                else :
                    common_edges_g1[edge1] = [ edge2 ]
                if common_edges_g2.has_key(edge2) : 
                    common_edges_g2[edge2].append(edge1)
                else :
                    common_edges_g2[edge2] = [ edge1 ]
    return common_edges_g1, common_edges_g2


def graphlist2dict(graph):
    """Convert graph:
    from: list of edges
    to:   dictionary of nodes"""
    edges_dict = {}
    for edge in graph:
        if not edges_dict.has_key(edge.v1):
            edges_dict[edge.v1] = set((edge.v2,))
        else:
            edges_dict[edge.v1].add(edge.v2)

        if not edges_dict.has_key(edge.v2):
            edges_dict[edge.v2] = set((edge.v1,))
        else:
            edges_dict[edge.v2].add(edge.v1)
    return edges_dict


def recursive_rm_node(edges_dict, node):
    if len(edges_dict[node]) >= 2:
        return
    if len(edges_dict[node]) == 0:
        edges_dict.pop(node)
        return
    nodep = list(edges_dict[node])[0]
    edges_dict.pop(node)
    edges_dict[nodep].remove(node)
    recursive_rm_node(edges_dict, nodep)


def rm_monogamous(edges_dict):
    """filter nodes not in a triangular pattern"""
    nodes = edges_dict.keys()
    for n in nodes:
        if edges_dict.has_key(n) and len(edges_dict[n]) < 2:
            recursive_rm_node(edges_dict,n)


def bfs2(edges_dict, queue=None, visited=None):
    """Breadth First Search iterator"""
    if visited is None:
        visited = set()
    
    if queue is None:
        queue = [edges_dict.keys()[0]]
    elif not queue:
        raise StopIteration
    
    # visit neighbouring nodes
    for n in edges_dict[queue.pop(0)]:
        if n not in visited:
            queue.append(n)
            visited.add(n)
            yield n

    # re-run bfs on the updated 'queue' and 'visited' lists 
    for n in bfs2(edges_dict, queue, visited):
            yield n
    

def allsubgraphs(edges_dict):
    subgraphs = []
    nodes = set(edges_dict.keys())        
    while len(nodes):
        subgraph = []
        start = list(nodes)[0]
        for node in bfs2(edges_dict, queue=[start]):    
            for connected_node in edges_dict[node]:
                subgraph.append(Edge(node, connected_node))
            nodes.remove(node)
        subgraphs.append(subgraph)
    return subgraphs

 
def build_subgraphs(common_edges):
    edges_dict = graphlist2dict(common_edges)
    rm_monogamous(edges_dict)
    nodes = set(edges_dict.keys())
    subgraphs = allsubgraphs(edges_dict)
    return allsubgraphs(edges_dict) 


def write_graph(graph, filename, append=False):
    """save graph in text file"""
    mode = 'a' if append else 'w'
    with open(filename, mode) as OUT:
        OUT.write("Node1\tx1\ty1\tz1\tNode2\tx2\ty2\tz2\n")
        for edge in graph:
            OUT.write(edge.tabrepr)
            

def build_graph_sub(pdb_fname, chain = None): 
    structure = load_PDB(pdb_fname)
    accessible_aa = select_aa_acessible(pdb_fname, chain)
    vertices = select_vertices(structure, accessible_aa)
    graph = build_graph(vertices)
    return graph


if __name__ == '__main__':
    graph1 = build_graph_sub('2CPK.pdb', chain='E')
    graph2 = build_graph_sub('3LCK.pdb', chain='A')
    common_edges1, common_edges2 = edges_comp(graph1, graph2)
    sgs1 = build_subgraphs(list(set(common_edges1)))
    sgs2 = build_subgraphs(list(set(common_edges2)))
    
    for i, sg in enumerate(sgs1, start=1):
        write_graph(sg, "2CPK_graphs_%02i.tsv" % i)
    for i, sg in enumerate(sgs2, start=1):
        write_graph(sg, "3LCK_graphs_%02i.tsv" % i)

