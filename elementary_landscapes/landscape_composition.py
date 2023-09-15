# Author: Mahakaran Sandhu 

import itertools
from landscape_model import make_NK
import networkx as nx
import numpy as np
import random
from scipy import sparse

# landscape generation

class NK_landscape: 
    """Class for NK landscapes"""
    def __init__(self, N, K, alphabet, distribution):
        self.landscape = make_NK(N, K, alphabet, distribution)
        self.N = N
        self.K = K
        self.alphabet = alphabet
        self.distribution = distribution        
        self.sequences = self.landscape[0]
        self.fitnesses = self.landscape[1]
        self.landscape_dict = {x:y for x,y in zip(self.sequences.flatten(), self.fitnesses.flatten())}

def all_genotypes(N, AAs):
    """Fills the sequence space with all possible genotypes."""
    return [''.join(i) for i in np.array(list(itertools.product(AAs, repeat=N)))]

# graph functions

def hamming_circle(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.
    
    """
    for positions in itertools.combinations(range(len(s)), n):
        for replacements in itertools.product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)

def graph_matrices(sequence_space, alphabet, incomplete=False): 
    """Returns adjacency, Laplacian and sequence space for a Hamming graph for a length N protein 
       over a specified alphabet."""
    
    seq_space = [''.join(list(i)) for i in sequence_space]
    seq_len    = len(seq_space[0])
    alph_len   = len(alphabet)
    
    
    lone = 0
    
    members    = set(seq_space)
    nodes      = {x:y for x,y in zip(seq_space, range(len(seq_space)))}
    adjacency  = sparse.lil_matrix((len(seq_space), len(seq_space)), dtype='int8') 
    
    for ind in range(len(sequence_space)):        
        seq = sequence_space[ind]     
        
        g=0
        
        for neighbor in hamming_circle(seq, 1,alphabet):
            if incomplete:
                try:
                    adjacency[ind,nodes[neighbor]]=1
                    g+=1
                except: 
                    pass
            else:                 
            
                adjacency[ind,nodes[neighbor]]=1
        if g==0:
            lone+=1
            
            
    degree    = (seq_len*(alph_len-1))*sparse.eye(len(seq_space))
    
    laplacian = degree - adjacency
    
    return laplacian.tocsc(), adjacency.tocsc(), degree.tocsc()

def cyclic_group_im(n): 
    """Returns cyclic group of order n in the complex numbers"""
    return np.array([np.exp((2*np.pi*1j*(i))/n) for i in range(n)])

    
def draw_radial_hierarchical(G,radius,alphabet, root=None): 
    """Draw a graph with hierarchical radial layout. Ensure that nodes have associated sequences in the attributes"""
    if root==None: 
        root = random.sample(graph.nodes(), k=1)[0]
    
    seed_seq       = G.nodes[root]['seq']
    max_dist       = len(seed_seq)
    neighborhoods  = [list(hamming_circle(seed_seq, n+1, alphabet=alphabet)) for n in range(max_dist)]
    positions_im   = [cyclic_group_im(len(n)) for n in neighborhoods]
    positions_1    = []
    r=radius
    for i in positions_im:
        positions_1.append(i*r)
        r+=radius
    positions_2    = [ list(zip(np.real(t), np.imag(t) )) for t in positions_1]
    
    #assign positions 
    positions_dict = [{x:y for x,y in zip(neighborhoods[i], positions_2[i])} for i in range(len(neighborhoods))]
    merged = {}
    for d in positions_dict:
        merged.update(d)
    for ind, i in enumerate(G.nodes()):
        try:
            G.nodes[i]['pos'] = merged[G.nodes[i]['seq']]
        except:
            pass
    G.nodes[root]['pos'] = (0,0)
    
    return G, nx.get_node_attributes(G,'pos'), merged


def write_to_csv(dictionary, filepath): 
    with open(filepath, 'w') as file:
        keys = list(dictionary.keys())
        vals = list(dictionary.values())
        for i in range(len(dictionary)): 
            file.write('{},{}\n'.format(keys[i], vals[i]))

if __name__ == "__main__":

    # initialise sequence space
    aas = "ACD"
    n = 5
    sequence_space = all_genotypes(n, aas)

    # get adjacency matrix and decomposition
    _, adj_mat, _ = graph_matrices(sequence_space, aas)
    eval, evec = np.linalg.eigh(adj_mat.toarray())

    # plot elementary landscapes
    evec_flipped = np.flip(evec, axis=1)
    graph = nx.from_numpy_matrix(adj_mat.toarray())
    for idx, i in enumerate(graph.nodes()):
        graph.nodes[i]['seq'] = sequence_space[idx]

    nx.write_edgelist(graph, "edge_list.csv")
    nx.write_gml(graph, "seq_graph.gml")

    # elementary landscapes
    # order 2
    order2 = evec_flipped[:, 10:50]
    # order 3
    order3 = evec_flipped[:, 50:130]

    # make composite
    order2_col_idx = range(order2.shape[1])
    order3_col_idx = range(order3.shape[1])
    order_2_3_combination_idx = list(itertools.product(
        order2_col_idx, order3_col_idx
    ))

    order2_landscape_ls = [{x:np.round(y, decimals=10) for x,y in zip(sequence_space, evec_flipped[:,i])} for i in range(10, 50)]
    order3_landscape_ls = [{x:np.round(y, decimals=10) for x,y in zip(sequence_space, evec_flipped[:,i])} for i in range(50, 130)]
    order2_3_comp_landscape_ls = [order2[:,i[0]] + order3[:,i[1]] for i in order_2_3_combination_idx]

    # example landscape and corresponding composite
    ## noting that this is an example of indices found in ipynb. 9 and 32 are the
    ## indices that combine to make the i=752 composite. Reproduction should be done 
    ## with the above code and searches in ipynb on ruggedness criteria of interest.
    ## (in our case Dirichlet energy)
    order2_landscape = order2_landscape_ls[9]
    order3_landscape = order3_landscape_ls[32]
    order2_3_comp_landscape = order2_3_comp_landscape_ls[752]

    write_to_csv(order2_landscape, 'order_2_landscape_node_fitnesses.csv')
    write_to_csv(order3_landscape, 'order_3_landscape_node_fitnesses.csv')
    write_to_csv(order2_3_comp_landscape, 'composite_landscape_node_fitnesses.csv')