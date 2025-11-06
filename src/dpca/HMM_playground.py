import numpy as np
from hmmlearn import hmm
import platform
from tqdm import tqdm
import gzip
from statistics import median
import lzma
import argparse
import shutil
from pathlib import Path
import random
import matplotlib.pyplot as plt
import _params
from _aux_functions import generate_sh_param_file, compress_file, build_maple_entry, build_maple_file

parser = argparse.ArgumentParser(description='Detect potential contaminated areas depending on several thresholds applied on the Viridian samples.')

parser.add_argument('--batch_id', type=int,
                    help='Batch id to process.')

args = parser.parse_args()

# Get the arguments
batch_id = args.batch_id

n_batch = _params.n_batch
het_thr = _params.het_thr
depth_thr = _params.depth_thr
prop_under_depth_thr = _params.prop_under_depth_thr
typical_depth_thr = _params.typical_depth_thr
n_masked_thr = _params.n_masked_thr
path_ref_seq = _params.path_ref_seq
inputTree = _params.inputTree
inputRates = _params.inputRates
n_diff_mut = _params.n_diff_mut
masked_max_dist = _params.masked_max_dist
masking_ratio = _params.masking_ratio
max_n_het_sites = _params.max_n_het_sites
param_term = _params.param_term

def generate_batchs_pairs(n_batchs):
    """Generate a 2D list with structure such that L[i][j] = sample j of batch i.
    Samples are pairs (path, read_name) where path is the path to the qc file.

    Args:
        n_batchs (int): Number of batches to generate.
    """
    path_vdn = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"
    # Generate a list of ids to chose samples
    batchs_samples_list = [[] for _ in range(n_batchs)]

    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        for i, line in enumerate(f):
            if i == 0:
                continue
            # Get the read_name and path
            read_name, path = line.strip().split("\t")
            path = path + "/qc.tsv.gz"
            batchs_samples_list[i % n_batchs].append((path, read_name))
            
            
        # for i, batch in tqdm(enumerate(batchs_samples_list), desc="Generating batches", total=n_batchs):
        #     with lzma.open(f"{data_storing_dir}/batches/batch_{i}.pkl.xz", "wb") as f:
        #         pickle.dump(batch, f)
        # Not enough memory to store all batches
    
    return batchs_samples_list


if __name__ == "__main__":
    
    print(platform.python_implementation(), platform.python_version())

    batchs_samples_list = generate_batchs_pairs(n_batch)
    batch_id = random.randint(0, n_batch - 1)
    sample_list = batchs_samples_list[batch_id]
    n_samples = len(sample_list)

    # Try a GAUSSIAN HMM
    # 4 Hidden states model accounting for the presence of contamination:
    # - State 0: Clean region (NN) ~ Coverage between 1000 and 2000
    # - State 1: Dropout of the main genome, contaminated region (DN) ~ Coverage between 20 and 150
    # - State 2: Dropout of main genome and contaminant (DD) ~ Coverage between 0 and 20
    # - State 3: Dropout of contaminant only (ND) ~ Coverage between 150 and 1000
    
    # Test with means: [1050, 85, 10, 350] and variances: [600, 65, 10, 200]
    # We set the means and variances manually to avoid bad initialization
    
    model_4 = hmm.GaussianHMM(
        n_components=4,
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
        n_iter=10,
        tol=1e-2,
        verbose=True,
        params="st",
        init_params="st",
    )
        
    model_4.means_ = np.array([[1050.0], [85.0], [10.0], [350.0]])
    model_4.covars_ = np.array([[600.0], [65.0], [10.0], [200.0]])
    
    # SECOND IDEA
    # 3 Hidden states model accounting for the presence of contamination:
    # Model based on sequences of cons_depth/clean_depth (heterozygosity)
    # - State 0: Clean region (NN) ~ 100% of reads support the consensus
    # - State 1: High heterozygous region (legitimate heterozygosity) ~ Mean around 50% of reads support the consensus
    # - State 2: Low heterozygosity region ~ Mean around 85% (~15% of reads come from the contaminant)
    
    # Test with means: [1, 0.5, 0.85] and variances: [0.05, 0.3, 0.1]
    
    model_3 = hmm.GaussianHMM(
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
        n_iter=10,
        tol=1e-2,
        verbose=True,
        params="st",
        init_params="st",
    )

    model_3.means_ = np.array([[1.0], [0.5], [0.85]])
    model_3.covars_ = np.array([[0.02], [0.3], [0.13]])
    
    # THIRD IDEA: In the 4 states model, the states NN and ND (0 and 3) are redundant for the purpose of contamination detection.
    # We could use a 3 states model where:
    # - State 0: Clean region (NN) ~ Coverage between 200 and 2000
    # - State 1: Full dropout (DD) ~ Coverage between 0 and 20
    # - State 2: Contaminated region (DN) ~ Coverage between 20 and 200
    # This would simplify the model and make it more robust to noise.
    
    # Test with means: [1000, 10, 110] and variances: [800, 10, 90]
    
    model_3_cont = hmm.GaussianHMM(
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
        n_iter=10,
        tol=1e-2,
        verbose=True,
        params="st",
        init_params="st",
    )
    
    model_3_cont.means_ = np.array([[1205.0], [951.0], [106.0]])
    model_3_cont.covars_ = np.array([[np.sqrt(305459.0)], [np.sqrt(838.0)], [np.sqrt(10404.0)]])
    
    # 3 states model is aimed at detecting on its own the contaminated sample, but does not detect the contaminated regions.
    # 4/3 states model is aimed at detecting the contaminated regions, and could be coupled with masking etc.

    
    for path, read_name in tqdm(sample_list,
                                desc="Processing samples",
                                total=n_samples,
                                unit="sample",
                                mininterval=300):
        
        clean_depth_list = []
        cons_depth_list = []
        heterozygosity_list = []
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
                    heterozygosity = cons_depth / clean_depth if clean_depth > 0 else 1
                    # Lowercase the reference and unmasked sequences
                    cons_nt = line[4].lower() # Consensus is Viridian's consensus
                except:
                    ref_pos = int(line[0])
                    clean_depth = 0
                    cons_depth = 0
                    cons_nt = line[4].lower()
                    heterozygosity = 1
                    # Skip the first and last lines where depths col are filled with .
                
                if (len(ref_pos_list) > 0) and (ref_pos == ref_pos_list[-1]):
                    # Insertion
                    unma_seq[-1] += cons_nt
                
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    clean_depth_list.append(clean_depth)
                    cons_depth_list.append(cons_depth)
                    heterozygosity_list.append(heterozygosity)
        
        # Apply the HMM to the clean_depth_list
        X_depth = np.array(clean_depth_list).reshape(-1, 1)
        X_het = np.array(heterozygosity_list).reshape(-1, 1)
        # Fit the model to the data
        try:
            model_4.fit(X_depth)
            hidden_states_4 = model_4.predict(X_depth)
            
            model_3_cont.fit(X_depth)
            hidden_states_3_cont = model_3_cont.predict(X_depth)
            
            model_3.fit(X_het)
            hidden_states_3 = model_3.predict(X_het)
        except:
            print(f"Error fitting HMM for sample {read_name}, skipping...")
            continue
        
        # Get the positions to mask (State 1 and State 2)
        positions_to_mask_4 = [ref_pos_list[i] for i in range(len(hidden_states_4)) if hidden_states_4[i] in [1, 2]]
        positions_to_mask_3 = [ref_pos_list[i] for i in range(len(hidden_states_3)) if hidden_states_3[i] in [1, 2]]
        
        # Print the parameters of the model after fitting
        print(f"Sample: {read_name}")
        print("Means:", model_4.means_.flatten())
        print("Variances:", model_4.covars_.flatten())
        print("Transition matrix:\n", model_4.transmat_)
        print("Start probabilities:", model_4.startprob_.flatten())
        print(f"Number of positions to mask: {len(positions_to_mask_4)} for model 4-states, {len(positions_to_mask_3)} for model 3-states (depth)")
        
        print("State sequence:")
        # Plot a figure of the state sequence
        plt.figure(figsize=(15, 5))
        plt.plot(clean_depth_list, label="Clean depth", color="blue", alpha=0.5)
        plt.scatter(range(len(hidden_states_4)), clean_depth_list, c=hidden_states_4, cmap="viridis", label="Hidden states")
        plt.xlabel("Position index")
        plt.ylabel("Clean depth")
        plt.title(f"Sample {read_name} - Clean depth and Hidden states")
        plt.colorbar(label="Hidden state")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"/nfs/research/goldman/anoufa/figs/hmm/hmm_depth_seq_{read_name}.png")
        
        # Same fig for 3-states alt model
        plt.figure(figsize=(15, 5))
        plt.plot(clean_depth_list, label="Clean depth", color="blue", alpha=0.5)
        plt.scatter(range(len(hidden_states_3_cont)), clean_depth_list, c=hidden_states_3_cont, cmap="viridis", label="Hidden states")
        plt.xlabel("Position index")
        plt.ylabel("Clean depth")
        plt.title(f"Sample {read_name} - Clean depth and Hidden states (3 states)")
        plt.colorbar(label="Hidden state")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"/nfs/research/goldman/anoufa/figs/hmm/hmm_depth_seq_{read_name}_3states.png")

        
        # Analysis of model 3
        print("Means (model 3):", model_3.means_.flatten())
        print("Variances (model 3):", model_3.covars_.flatten())
        print("Transition matrix (model 3):\n", model_3.transmat_)
        print("Start probabilities (model 3):", model_3.startprob_.flatten())
        # Plot a figure of the state sequence
        plt.figure(figsize=(15, 5))
        plt.plot(heterozygosity_list, label="Heterozygosity (cons_depth/clean_depth)", color="orange", alpha=0.5)
        plt.scatter(range(len(hidden_states_3)), heterozygosity_list, c=hidden_states_3, cmap="viridis", label="Hidden states (model 3)")
        plt.xlabel("Position index")
        plt.ylabel("Heterozygosity")
        plt.title(f"Sample {read_name} - Heterozygosity and Hidden states (model 3)")
        plt.colorbar(label="Hidden state")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"/nfs/research/goldman/anoufa/figs/hmm/hmm_het_seq_{read_name}.png")
        
        # Print the number of positions in each state for both models
        for i in range(model_4.n_components):
            print(f"Model 4 - State {i}: {np.sum(hidden_states_4 == i)} positions")
        for i in range(model_3.n_components):
            print(f"Model 3 - State {i}: {np.sum(hidden_states_3 == i)} positions")
        for i in range(model_3_cont.n_components):
            print(f"Model 3 cont - State {i}: {np.sum(hidden_states_3_cont == i)} positions")
            
        print("Param of 3 states contamination model:")
        print("Means (model 3 cont):", model_3_cont.means_.flatten())
        print("Variances (model 3 cont):", model_3_cont.covars_.flatten())
        print("Transition matrix (model 3 cont):\n", model_3_cont.transmat_)
        print("Start probabilities (model 3 cont):", model_3_cont.startprob_.flatten())
        
        break
        
                
    
    
                            