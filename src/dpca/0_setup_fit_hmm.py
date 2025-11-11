from tqdm import tqdm
import gzip
import lzma
import argparse
import shutil
from pathlib import Path
import random
from hmmlearn import hmm
import numpy as np
import _params
from _aux_functions import generate_sh_param_file, compress_file, build_maple_entry, build_maple_file, update_params_file

startprob = _params._startprob
transmat = _params._transmat
means = _params._means
covars = _params._covars
samples_dir = _params.samples_dir

def generate_sample_data(n_samples):
    """Generate a 2D list with structure such that L[i][j] = sample j of batch i.
    Samples are pairs (path, read_name) where path is the path to the qc file.

    Args:
        n_batchs (int): Number of batches to generate.
    """
    path_vdn = samples_dir
    # Generate a list of ids to chose samples
    batchs_samples_list = []
    
    # DRAW n_samples integer between 0 and 4,900,000
    drawn_samples = random.sample(range(4900000), n_samples)

    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        for i, line in enumerate(f):
            if i == 0:
                continue
            
            if i in drawn_samples:
                # Get the read_name and path
                read_name, path = line.strip().split("\t")
                path = path + "/qc.tsv.gz"
                batchs_samples_list.append((path, read_name))
    
    return batchs_samples_list

def setup_fit_hmm(sample_data, lengths):
    """
    Given a sample data dictionary, set up and fit an HMM to the depth data.
    """
    # THIRD IDEA: In the 4 states model, the states NN and ND (0 and 3) are redundant for the purpose of contamination detection.
    # We could use a 3 states model where:
    # - State 0: Clean region (NN) ~ Coverage between 200 and 2000
    # - State 1: Full dropout (DD) ~ Coverage between 0 and 20
    # - State 2: Contaminated region (DN) ~ Coverage between 20 and 200
    # This would simplify the model and make it more robust to noise.
    
    # Test with means: [1000, 10, 110], variances: [800, 10, 90] and startprob [0.0, 1.0, 0.0]

    hmm_model = hmm.GaussianHMM(
        n_components=3,
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
        n_iter=20,
        tol=1e-2,
        verbose=False,
        params="t",
        init_params="t",
    )
    
    hmm_model.means_ = np.array(means).reshape(-1, 1)
    hmm_model.covars_ = np.array(covars).reshape(-1, 1)
    hmm_model.startprob_ = np.array(startprob, dtype=float).ravel()
    
    # Fit the model to the data
    X = np.array(clean_depth_list).reshape(-1, 1)
    hmm_model.fit(X, lengths=lengths)
    
    print("Model fitted.")
    
    # SAVE THE MODEL PARAMETERS IN PARAMS.PY
    
    model_params = {
        "_startprob": hmm_model.startprob_.tolist(),
        "_transmat": hmm_model.transmat_.tolist(),
        "_means": hmm_model.means_.flatten().tolist(),
        "_covars": hmm_model.covars_.flatten().tolist(),
    }

    params_path = "/nfs/research/goldman/anoufa/src/dpca/_params.py"
    print(f"Parameters to be saved: {model_params}")
    update_params_file(params_path, model_params)


if __name__ == "__main__":
    n_samples = 5000
    sample_tuple_list = generate_sample_data(n_samples)
    # Convert the list of tuples to a list of sequences
    clean_depth_list = []
    for path, read_name in tqdm(sample_tuple_list,
                                desc="Processing samples",
                                total=n_samples,
                                unit="sample",
                                mininterval=300):
        cons_depth_list = []
        ref_pos_list = []
        unma_seq = []
        
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
                    # Lowercase the reference and unmasked sequences
                    cons_nt = line[4].lower() # Consensus is Viridian's consensus
                except:
                    ref_pos = int(line[0])
                    clean_depth = 0
                    cons_depth = 0
                    cons_nt = line[4].lower()
                    # Skip the first and last lines where depths col are filled with .
                
                if (len(ref_pos_list) > 0) and (ref_pos == ref_pos_list[-1]):
                    # Insertion
                    unma_seq[-1] += cons_nt
                
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    clean_depth_list.append(clean_depth)
                    cons_depth_list.append(cons_depth)
    
    print("List of clean depths generated, now fitting the HMM model...", flush=True)
    lengths = [29903]*n_samples

    setup_fit_hmm(clean_depth_list, lengths)