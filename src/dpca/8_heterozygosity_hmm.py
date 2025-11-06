# This script should go over the tsv results of MAPLE sample placement and
# generate the distributions of the branch lengths for the three types of sequences
# (viridian, masked, random).

# The goal is to see whether the masked samples are closer to the tree than the viridian ones or not.

import pickle
from pathlib import Path
import gzip
from tqdm import tqdm
from collections import Counter
import _params
import argparse
import csv
import random
from _aux_functions import update_params_file, generate_sample_list
from multiprocessing import Pool
from os import cpu_count
import os
import numpy as np
import pandas as pd
import sys
sys.path.append("/nfs/research/goldman/anoufa/src")
from dpca.hmmlearn.src.hmmlearn import hmm
import lzma


parser = argparse.ArgumentParser(
    description="Process the concatenated output of the MAPLE sample placement."
)

parser.add_argument(
    "--masked_or_random",
    type=str,
    help="Either masked or random, to process both in parallel",
)

args = parser.parse_args()

# Get the arguments
masked_or_random = args.masked_or_random
path_psmp_masked = _params.masked_with_scores
path_psmp_random = _params.random_with_scores
param_term = _params.param_term
num_cores = _params.num_cores
startprob = _params._startprob
# transmat = _params._transmat
# means = _params._means
covars = _params._covars




def init_hmm(n_comp, covars, startprob):
    # Initialize the HMM model with the parameters from params.py
    hmm_model = hmm.GaussianHMM(
    n_components=n_comp,
    covariance_type="diag",
    min_covar=1e-3,
    startprob_prior=1.0,
    transmat_prior=1.0,
    means_prior=0,
    means_weight=0,
    covars_prior=1e-2,
    covars_weight=1,
    algorithm="viterbi",
    random_state=7,
    n_iter=50,
    tol=1e-2,
    verbose=False,
    params="tm",
    init_params="",
    )

    covars = covars[:n_comp]
    startprob = startprob[:n_comp]
    # Initialize the transmat with uniform probabilities
    transmat = np.full((n_comp, n_comp), 1.0 / n_comp)
    means = [0.0 for _ in range(n_comp)]
    
    hmm_model.means_ = np.array(means).reshape(-1, 1)
    hmm_model.covars_ = np.array(covars).reshape(-1, 1)
    hmm_model.startprob_ = np.array(startprob, dtype=float).ravel()
    hmm_model.transmat_  = np.array(transmat, dtype=float).reshape(n_comp, n_comp)

    return hmm_model

if __name__ == "__main__":

    # Check that we are correctly using the local hmmlearn
    print(f"Using hmmlearn from: {hmm.__file__}")
    if masked_or_random == "masked":
        path_store = path_psmp_masked

    elif masked_or_random == "random":
        path_store = path_psmp_random

    df = pd.read_csv(path_store, sep="\t")
    sample_names = set(df["sample_name"].tolist())
    print(f"Generating sample list: {len(sample_names)} {masked_or_random} samples to process")
    sample_list = generate_sample_list(sample_names)
    print(f"Sample list generated: {len(sample_list)} samples to process")
    hmm_model_2 = init_hmm(n_comp=2, covars=covars, startprob=startprob)
    hmm_model_3 = init_hmm(n_comp=3, covars=covars, startprob=startprob)
    is_cont_dict = {}
    
    # Subsample sample list
    # sample_list = random.sample(sample_list, 200)
    
    for path, read_name in tqdm(sample_list,
                                desc="Processing samples",
                                unit="sample",
                                mininterval=120):
        clean_depth_list = []
        cons_depth_list = []
        ref_pos_list = []
        unma_seq = []
        heterozygosity_list = []
        minor_alleles_list = []
        
        # unzip and open the file
        with gzip.open(path, "rb") as f:
            # Read the lines and create the needed values to apply the thresholds
            # List of all the clean_depth, list of het_prop, median_depth
            
            # b'Ref_pos\tRef_nt\tCons_pos\tCons_nt\tMasked_cons_nt\tAmplicon\tPrimer\tMask\tTotal_depth\tClean_depth\tCons_depth\tA\ta\tC\tc\tG\tg\tT\tt\tI\ti\tD\td\tX_A\tX_a\tX_C\tX_c\tX_G\tX_g\tX_T\tX_t\tX_I\tX_i\tX_D\tX_d\n'
            
            # Iterate over the lines in the file
            for line in f:
                # ignore header line
                if line.startswith(b"Ref_pos"):
                    continue
                # Decode the line and split it into columns
                line = line.decode("utf-8").strip().split("\t")
                # Extract the values from the columns
                try:
                    # Extract the values from the columns
                    ref_pos = int(line[0])
                    clean_depth = int(line[9])
                    cons_depth = int(line[10])
                    # Find the second highest nucleotide count
                    nt_counts_dict = {"A": int(line[11]) + int(line[12]), "C": int(line[13]) + int(line[14]), "G": int(line[15]) + int(line[16]), "T": int(line[17]) + int(line[18])}
                    # Set the max to 0
                    nt_counts_dict[max(nt_counts_dict, key=nt_counts_dict.get)] = 0
                    minor_allele = max(nt_counts_dict, key=nt_counts_dict.get)
                    het_depth = nt_counts_dict[minor_allele]/clean_depth if clean_depth > 0 else 0
                    # Lowercase the reference and unmasked sequences
                    cons_nt = line[4].lower() # Consensus is Viridian's consensus
                except:
                    ref_pos = int(line[0])
                    clean_depth = 0
                    cons_depth = 0
                    cons_nt = line[4].lower()
                    het_depth = 0
                    minor_allele = 'N'
                    
                    # Skip the first and last lines where depths col are filled with .
                
                if (len(ref_pos_list) > 0) and (ref_pos == ref_pos_list[-1]):
                    # Insertion
                    unma_seq[-1] += cons_nt
                
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    clean_depth_list.append(clean_depth)
                    cons_depth_list.append(cons_depth)
                    heterozygosity_list.append(het_depth)
                    minor_alleles_list.append(minor_allele)

        # Apply heterozygosity HMM to identify the contaminant percentage
        hmm_model_2.fit(np.array(heterozygosity_list).reshape(-1, 1))
        # Check the convergence history of the sample
        print(f"Means for sample {read_name} with 2-states HMM: {hmm_model_2.means_.flatten()}")
        if not hmm_model_2.monitor_.converged:
            # Try the 3-states HMM
            hmm_model_3.fit(np.array(heterozygosity_list).reshape(-1, 1))
            hmm_model = hmm_model_3
            print(f"Sample {read_name}: 2-states HMM did not converge, but 3-states HMM: {hmm_model_3.monitor_.converged}")
        else:
            hmm_model = hmm_model_2
        # Two solutions: either both states are ~0 and we say no contamination
        # Or one state is high het and we should store the positions and the nucleotides of the minor alleles to check if they match the contaminant
        tol = 0.05
        if all(mean < tol for mean in hmm_model.means_.flatten()):
            # No contamination detected
            is_cont_dict[read_name] = (False, {})
        
        else:
            hidden_states = hmm_model.predict(np.array(heterozygosity_list).reshape(-1, 1))
            high_mean_state = np.argmax(hmm_model.means_)
            het_minor_alleles_dict = {}
            for i, state in enumerate(hidden_states):
                if state == high_mean_state:
                    het_minor_alleles_dict[ref_pos_list[i]] = minor_alleles_list[ref_pos_list[i]]
            # Store the positions and the minor alleles
            is_cont_dict[read_name] = (True, het_minor_alleles_dict)
            
    # Store the contamination dictionary
    path_is_cont_dict = Path(f'/nfs/research/goldman/anoufa/data/MAPLE_output/final_filter/is_cont_dict_{masked_or_random}.pickle')
    with open(path_is_cont_dict, 'wb') as handle:
        pickle.dump(is_cont_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Update the params file with the new path
    params_path = "/nfs/research/goldman/anoufa/src/dpca/_params.py"

    update_params_file(
        params_path,
        {
            f"is_cont_dict_{masked_or_random}": str(
                path_is_cont_dict
            )
        },
    )
    
    # Print the number of contaminated samples detected
    n_contaminated = sum(1 for v in is_cont_dict.values() if v[0])
    print(f"Number of contaminated samples detected in {masked_or_random}: {n_contaminated} out of {len(is_cont_dict)}")
    
    # Print an example entry that's true
    for sample, (is_cont, het_dict) in is_cont_dict.items():
        if is_cont:
            print(f"Example contaminated sample in {masked_or_random}: {sample} with {het_dict} heterozygous positions")
            break
        
    # Look for the best contaminant of the sample and check if the heterozygous positions match
