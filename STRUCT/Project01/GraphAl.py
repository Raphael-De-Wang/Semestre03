#!env python


#### Part 1 : construction of graph

def load_PDB():
    pass

def select_aa_acessible():
    pass


def select_vertices(atoms_mtx,typeChimi=["CB","CA"]):
    for a in atoms_mtx :
        if typeChi(a) in typeChimi :
            pass
    return vertices

def build_mtx_dist(vertices_list):
    N = len(vertices_list)
    for i in range(N):
        for j in range(i+1, N):
            compute dist
    return mtx

def weight(i,j):
    S1S2 = distance entre CA de i et CA de j
    S1E2
    S2E1
    q = angle entre i et j
    return (,,,)


def build_graph(vertices_list, Dmax=8): #Dmax angstrom
    mtx = build_mtx_dist(vertices_list)
    mtx[mtx<=Dmax]

    graph = []
    for i in range(N):
        for j in range(i+1, N):
            compute dist
            if dist <= Dmax:
                graph.append(i, j, weight(i,j))
                
    return graph     # [(node1,node2, weight<vector>) for i in N ]
    
def write_graph() : #
    pass


#### Part 2 : comparison of graph

def residu_ident(edge1,edge2):
    pass
def weight_ident(edge1,edge2):
    pass

def edges_comp(graph1, graph2):
    common_edges_g1 = []
    common_edges_g2 = []
    for edge1 in grapĥ1:
        for edge2 in graph2:
            if residu_ident(edge1,edge2) and weight_ident(edge1,edge2):
                common_edges_g1.append(edge1)
                common_edges_g2.append(edge2)
    return (common_edges_g1,common_edges_g2)

def breadth_first_algo(graph):
    for ...
        for ..
            yield element
            

def build_subgraph(graph,common_edges):
    # fiter nodes no clique >= 3
    for node in breadth_first_algo(graph, common_edges):
        pass
    return subgraph_list
    

