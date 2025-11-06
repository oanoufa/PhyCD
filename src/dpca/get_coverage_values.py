# Generate distributions over the whole dataset

import platform
from tqdm import tqdm
import gzip
from statistics import median
import lzma
import argparse
import shutil
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import random


def generate_sample_data(n_samples):
    """Generate a 2D list with structure such that L[i][j] = sample j of batch i.
    Samples are pairs (path, read_name) where path is the path to the qc file.

    Args:
        n_batchs (int): Number of batches to generate.
    """
    path_vdn = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"
    # Generate a list of ids to chose samples
    batchs_samples_list = []
    
    # DRAW n_samples integer between 0 and 4,900,000
    drawn_samples = random.sample(range(4900000), n_samples)

    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        for i, line in tqdm(enumerate(f), desc="Processing lines", unit="line",
                            mininterval=180):
            if i == 0:
                continue
            
            if i in drawn_samples:
                # Get the read_name and path
                read_name, path = line.strip().split("\t")
                path = path + "/qc.tsv.gz"
                batchs_samples_list.append((path, read_name))
    
    return batchs_samples_list

if __name__ == "__main__":
    # Generate the batchs
    print(platform.python_implementation(), platform.python_version())
    cov_val_path = Path("/nfs/research/goldman/anoufa/data/gen_metrics/coverage_values.pkl")

    if cov_val_path.exists():
        print(f"Coverage values file {cov_val_path} already exists. Exiting.")
        coverage_values = pickle.load(open(cov_val_path, "rb"))
        print(f"Loaded existing coverage values of length {len(coverage_values)}")
    else:
        print(f"Generating new coverage values and storing them at {cov_val_path}")
        n_samples = 10000
        sample_data = generate_sample_data(n_samples)
        coverage_values = []
        # Store the coverage values of all those samples in a list
        for path, read_name in tqdm(sample_data,
                                    desc="Processing samples",
                                    total=n_samples,
                                    unit="sample",
                                    mininterval=30):
            cons_depth_list = []
            ref_pos_list = []
            unma_seq = []
            clean_depth_list = []

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
        
            coverage_values.extend(clean_depth_list)
        # Store the coverage_values using pickle
        with open(cov_val_path, "wb") as f:
            pickle.dump(coverage_values, f)
        print(f"Stored coverage_values of {n_samples} samples at {cov_val_path}")
    
    print(f"Mean coverage: {sum(coverage_values)/len(coverage_values)}, Median coverage: {median(coverage_values)}")

    # Plot a histogram of the coverage values
    plt.figure(figsize=(10, 6))
    sns.histplot(coverage_values, bins=100, kde=False)
    plt.title("Histogram of coverage values")
    plt.xlabel("Coverage values")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.tight_layout()
    plt_path = Path("/nfs/research/goldman/anoufa/figs/coverage_histogram.png")
    plt.savefig(plt_path)
    print(f"Saved histogram plot at {plt_path}")
    
    # Fit a Gaussian Mixture Model to the coverage values
    from sklearn.mixture import GaussianMixture
    import numpy as np
    
    # THREE COMPONENTS GMM
    # Subsample the coverage values for fitting
    # Clip coverage values to max 1600
    coverage_values = [min(cv, 1600) for cv in coverage_values]
    sample_size = min(100000, len(coverage_values))
    coverage_sample = random.sample(coverage_values, sample_size)
    X = np.array(coverage_sample).reshape(-1, 1)
    gmm = GaussianMixture(n_components=3, random_state=42)
    gmm.fit(X)
    print("Fitted GMM to coverage values.")
    # PRINT THE MEANS AND COVARIANCES OF THE GMM COMPONENTS
    for i in range(gmm.n_components):
        print(f"Component {i}: Mean = {gmm.means_[i][0]}, Covariance = {gmm.covariances_[i][0][0]}, Weight = {gmm.weights_[i]}")
    
    # Plot the GMM components over the histogram
    x = np.linspace(0, max(coverage_values), 1000).reshape(-1, 1)
    logprob = gmm.score_samples(x)
    responsibilities = gmm.predict_proba(x)
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]
    plt.figure(figsize=(10, 6))
    sns.histplot(coverage_values, bins=100, kde=False, stat="density",
                    label="Coverage Histogram")
    plt.plot(x, pdf, '-k', label='GMM Total')
    plt.plot(x, pdf_individual, '--', label='GMM Components')
    plt.title("GMM Fit to Coverage Values")
    plt.xlabel("Coverage values")
    plt.ylabel("Density")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    gmm_plt_path = Path("/nfs/research/goldman/anoufa/figs/coverage_gmm_fit.png")
    plt.savefig(gmm_plt_path)
    print(f"Saved GMM fit plot at {gmm_plt_path}")