from argparse import ArgumentParser
import numpy as np
import pandas as pd
from pathlib import Path
import random
from scipy import spatial
from sklearn.neighbors import NearestNeighbors

def get_ohe(seq_ls):
    '''
    Get the one-hot embedding for a list of sequences.
    
    Assumes sequences are aligned and any non-canonical amino acids have been
    handled.
    '''
    # make ohe dict
    aa_ls = set([aa for seq in seq_ls for aa in list(seq)])
    n_aas = len(aa_ls)
    ohe_dict = {}
    for i, aa in enumerate(aa_ls):
        arr_idx = np.zeros((n_aas))
        arr_idx[i] = 1
        ohe_dict[aa] = arr_idx

    # convert seqs to ohe
    ohe_seqs = []
    for seq in seq_ls:
        seq_ohe = np.empty(())
        for aa in seq:
            seq_ohe = np.hstack((seq_ohe, ohe_dict[aa]))
        ohe_seqs.append(seq_ohe)

    return np.stack(ohe_seqs)

def get_dirichlet_energy(ohe_arr, y):
    '''
    Get the normalized Dirichlet energy given the OHE and signals (y).
    '''
    # Make kNN
    dist_arr = spatial.distance_matrix(ohe_arr, ohe_arr)
    k = int(np.sqrt(len(y)))
    knn_fn = NearestNeighbors(
        n_neighbors=k,
        metric="precomputed"
    ).fit(dist_arr)

    # Determine graph Laplacian
    adj_mat = knn_fn.kneighbors_graph(dist_arr).toarray()
    ## remove self connections
    adj_mat -= np.eye(adj_mat.shape[0])
    ## make adjacency matrix symmetric
    adj_mat = adj_mat + adj_mat.T
    adj_mat[adj_mat > 1] = 1
    diag_mat = np.diag(np.sum(adj_mat, axis=0))

    laplacian = diag_mat - adj_mat
    dir_en = (y @ laplacian) @ y.T
    norm_dir_en = dir_en/len(y)
    return norm_dir_en

def bootstrap_dirichlet_energy(ohe_arr, y, n_samples=1000, subsample=0.95):
    '''
    Calculate the Dirichlet energy over subgraphs. 
    '''
    dir_en_ls = []

    for i in range(n_samples):
        # make random subsample of graph
        rand_idx = random.sample(list(range(len(y))), int(len(y) * subsample))
        y_subsample = y[rand_idx]
        ohe_arr_subsample = ohe_arr[rand_idx, :]
        norm_dir_en = get_dirichlet_energy(ohe_arr_subsample, y_subsample)
        dir_en_ls.append(norm_dir_en)

    mean_dir_en = np.mean(dir_en_ls)
    stdev_dir_en = np.mean(dir_en_ls)

    return mean_dir_en, stdev_dir_en


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--training_data', required=True, type=Path)
    args = parser.parse_args()
    
    # df containing LacI sequences & log enrichment values
    df = pd.read_csv("./lacI_enrichment.csv") 
    
    log_enrich = df.Log_enrichment.to_numpy()
    exp_log_enrich = df.Expected_log_enrichment.to_numpy()

    seq_ls = df.Sequence.tolist()
    ohe_arr = get_ohe(seq_ls)

    # actual dirichlet energy
    mean_dir_en, std_dir_en = bootstrap_dirichlet_energy(
        ohe_arr, 
        log_enrich
    )
    print("Actual Dirichlet energy:")
    print(f"mean = {mean_dir_en}, stdev = {std_dir_en}")

    # expected dirichlet energy
    exp_mean_dir_en, exp_std_dir_en = bootstrap_dirichlet_energy(
        ohe_arr, 
        exp_log_enrich
    )
    print("Expected Dirichlet energy:")
    print(f"mean = {exp_mean_dir_en}, stdev = {exp_std_dir_en}")