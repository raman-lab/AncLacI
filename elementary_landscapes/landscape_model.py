# Author: Mahakaran Sandhu 

import itertools
import numpy as np

def genEpiNet(N, K): 
    """Generates a random epistatic network for a sequence of length
    N (int) with, on average, K connections (int)"""
    return {
        i: sorted(np.random.choice(
            [n for n in range(N) if n != i], 
            K, 
            replace=False
        ).tolist() + [i])
        for i in range(N)
    }

def fitness_i(sequence, i, epi, mem, distribution, **kwargs):
    """Assigns a random fitness value to the ith amino acid that interacts with K other positions in a sequence.
    Fitnesses are drawn from a specified distribution (e.g. numpy.random.normal), with **kwargs passed to the 
    distribution sampling function. """

    key_0 = tuple(zip(epi[i], sequence[epi[i]]))    #we use the epistasis network to work out what the relation is  
    key   = (i, 0) + key_0

    if key not in mem:
        mem[key] = distribution(**kwargs)           #then, we assign a random number from distribution to this interaction
    return mem[key]

def fitness(sequence, epi, mem, distribution, residue_fitnesses=False, **kwargs):
    """Obtains a fitness value for the entire sequence by obtaining the mean over individual amino acids contributions."""
    per_residue = [fitness_i(sequence, i, epi, mem, distribution, **kwargs) for i in range(len(sequence))]
    if residue_fitnesses: 
        return np.mean(per_residue), per_residue
    else: 
        return np.mean(per_residue)

def min_max_scaler(x): 
    """Scales data in 1D array x to the interval [0,1]."""
    return (x-min(x))/(max(x)-min(x))

def make_NK(N, K, AAs, distribution, epi_net=None, minmax=True, residue_fitnesses=False, **kwargs): 
    """Make NK landscape with above parameters"""
    
    assert N>K, "K is greater than or equal to N. Please ensure K is strictly less than N."
    f_mem             = {}
    if epi_net is not None: 
        epistasis_network = epi_net
    else: 
        epistasis_network = genEpiNet(N, K)
    seq_space         = np.array(list(itertools.product(AAs, repeat=N)))

    if residue_fitnesses:
        fitness_tuple     = np.array([fitness(i, epi=epistasis_network, mem=f_mem, distribution=distribution, residue_fitnesses=residue_fitnesses, **kwargs) for i in seq_space])
        fitnesses         = np.array([float(i) for i in fitness_tuple[:,0]])
        fitnesses_i       = fitness_tuple[:,1]

    else: 
        fitnesses         = np.array([fitness(i, epi=epistasis_network, mem=f_mem, distribution=distribution, **kwargs) for i in seq_space])
        

    if minmax: 
        fitnesses     = min_max_scaler(fitnesses)


    seq_space         = np.array([''.join(list(i)) for i in itertools.product(AAs, repeat=N)]) #recalculate seq_space so its easier to concat
    
    fitnesses         = fitnesses.reshape((len(fitnesses),1))
    seq_space         = seq_space.reshape((len(seq_space),1))

    landscape         = (seq_space, fitnesses)
    
    if residue_fitnesses: 
        return landscape, epistasis_network, fitnesses_i, f_mem
    else: 
        return landscape