#################################################################################################################
# Mixed infection estimator                                                                                     #
# An algorithm for detection of mixed infections in bacterial whole genome sequence data.                       #
# The algorithm analyses a set of defined variable sites for evidence of mixed infection with up to 2 strains,  #
# if a mixed infection is present, the relative proportions of the dominant and minor haplotypes are estimated. #
# A database of known haplotypes is used to identify the most likely dominant and minor haplotype.              #
#                                                                                                               #
# david.eyre@ndm.ox.ac.uk                                                                                       #
# 04 March 2013                                                                                                 #
#                                                                                                               #
# This file is a python version of the original R code.                                                         #
#################################################################################################################

### LIBRARIES ###
import math
import numpy as np
from scipy.optimize import minimize
from scipy.special import logit
from itertools import product
import random
import os
from tqdm import tqdm
from pathlib import Path
from time import time
import argparse
from _aux_functions import update_params_file, load_pickle_dict, smart_open, compress_file
import plotly.figure_factory as ff

### PARSE ARGUMENTS ###
parser = argparse.ArgumentParser(
    description="Apply Eyre model to the potentially contaminated samples."
)
parser.add_argument("--masked_or_random", type=str,
                    help="Either masked or random, to process both in parallel")
parser.add_argument('--data_dir', type=str,
                    help='Directory containing all the files generated during the pipeline.')
parser.add_argument('--param_term', type=str,
                    help='str built from the initial parameters used to recognize the files.')
parser.add_argument('--compress', type=int,
                    help='1 or 0, enables compression.')
args = parser.parse_args()

param_term = args.param_term
data_dir = args.data_dir
masked_or_random = args.masked_or_random
compress = args.compress

### OPTIONS ###

#path to directory with base counts
path = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/")
# path = Path("/nfs/research/goldman/anoufa/src/test_eyre_model/Dataset_S1/")

#file with database of known haplotypes - tab delimited format with 2 columns: id, sequence
#sequences supplied as single string, each must be aligned and of the same length with no missing data
sequences_file = path / "haplotype_sequence.txt"

#file with population frequency of each haplotype
#tab delimited with 2 columns: id, frequency
sequences_freq_file = path / "haplotype_frequency.txt"

#file with list of file names containing the base counts with new line for each file
path_master_file = path / "filenames_list.txt"
with open(path_master_file)  as f:
    files = f.read().splitlines()
#files is now a list of file names

#file with pairs of haplotypes to consider for matching:
# dict with structure {
#     "sample_id1": {"contaminatnt_id1", "contaminatnt_id2", ...},
#     "sample_id2": {"contaminatnt_id3", "contaminatnt_id4", ...},
#     ...
# }
pairs_dict_file = path / f"sample_contaminant_pairs_{masked_or_random}.pickle"
variable_sites_per_sample_dict_path = path / f"variable_sites_per_sample_{masked_or_random}.pickle"

#file to record output
outlog = Path(f"{data_dir}/7/eyre_model_{masked_or_random}_output.tsv")
if not outlog.exists():
    outlog.parent.mkdir(parents=True, exist_ok=True)
    outlog.touch()
do_bootstrap = False  #whether to perform bootstrap to obtain confidence intervals
t0 = time()

# small value to use in place of 0 probabilities to avoid log(0)
err_correc = 1e-300

#number of cores to use for bootstrapping on Mac / Linux, ignored on Windows
bs_cores = 4

#number of bootstrap iterations to run for parameter estimates
bs_iter = 1000

#prior probability of a sample being contaminated
pi_con = 1e-3

### END OPTIONS ###

### CONSTANTS ###
DIPLOID = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
BASES = ["A", "C", "G", "T"]


### FUNCTIONS ###
#functions for defining m, mixture proportion and d, between haplotype diversity by maximum likelihood

#define functions to keep values between 0.5 and 1
def get_logit2(p):
    """Compute logit2 transformation.

    Args:
        p (float): Probability value between 0.5 and 1.

    Returns:
        float: Logit2 transformed value.
    """
    return math.log((2*p-1)/(1-(2*p-1)))
        

def get_inv_logit2(p):
    """Compute inverse logit2 transformation.

    Args:
        p (float): Logit transformed value.

    Returns:
        float: Probability value between 0 and 1.
    """
    return 0.5 + (math.exp(p) / (1 + math.exp(p))) / 2

#define functions to keep values between 0 and 1
def get_logit(p):
    """Compute logit transformation.

    Args:
        p (float): Probability value between 0 and 1.

    Returns:
        float: Logit transformed value.
    """
    return math.log(p/(1-p))

def get_inv_logit(p):
    """Compute inverse logit transformation.

    Args:
        p (float): Logit transformed value.

    Returns:
        float: Probability value between 0 and 1.
    """
    return math.exp(p) / (1 + math.exp(p))

#several functions to generate the likelihood of m, the mixture proportion
def get_p_b(A, m, b, e):
    """returns p(b|A,m,e) where A = {AA, AC, ..., TT} and b is a base in a single read, and e is the read error probability
    
    Args:
        A (str): Diploid genotype (e.g., "AA", "AC", ..., "TT").
        m (float): Mixture proportion.
        b (str): Base in a single read.
        e (float): Read error probability. (epsilon)
    """
    p_b_1 = m * e/3  # p where major ≠ b (so if b was read on major it's an error)
    p_b_2 = (1-m) * e/3  # p where minor ≠ b (so if b was read on minor it's an error)
    if b == A[0]:
        p_b_1 = m * (1-e)  # p where major = b (so if b was read on major there's no error)
    if b == A[1]:
        p_b_2 = (1-m) * (1-e)  # p where minor = b (so if b was read on minor there's no error)
    p_b = p_b_1 + p_b_2
    return p_b

def get_p_b_diploid(m, epsilon):
    """returns a matrix of p(b|A,m,e) for all values of b {A, C, G, T} - rows and values of A - columns
    """
    # e is set by constant above
    mat = np.empty((16, 4))
    for i, A in enumerate(DIPLOID):
        mat[i, :] = [get_p_b(A, m, b, epsilon) for b in BASES]
    return mat

def get_p_site(site, m, bc, mat):
    """returns p(B1 to Bn | A, m, e) where B1 to Bn are all reads at a given site
    - product of p(b|A,m,e) across all reads - taken from base counts
    where mat is a 16x4 matrix from get.p.b.diploid
    
    Output is a vector of length 16 with p(B1 to Bn | A_j, m, e) for each value of A_j
    """
    counts = bc[site, 1:5]  # ensure four BASES
    return np.exp(np.sum(np.log(mat) * counts, axis=1))

def get_p_site_for_A(site, m, bc, mat, A_j):
    """returns p(B1 to Bn | A_j, m, e) where B1 to Bn are all reads at a given site
    - product of p(b|A_j,m,e) across all reads - taken from base counts
    where mat is a 16x4 matrix from get.p.b.diploid
    and A_j is the specific value of A (0..15 index) to consider
    """
    return get_p_site(site, m, bc, mat)[A_j]

def logliki_m_epsilon(parm, bc, A_known):
    """
    Log-likelihood when haplotype pair μ is known.
    parm = [logit2(m), logit1(epsilon)]
    bc = base counts matrix
    A_known = array of known A values for each site (0..15 index)
    A_known is sequences_pairs_mat for the specific haplotype pair being considered
    """
    m = get_inv_logit2(logit_m)
    epsilon = get_inv_logit(logit_epsilon)
    
    mat = get_p_b_diploid(m, epsilon)

    log_lik = 0.0

    for site in range(bc.shape[0]):
        # compute P(B_site | this specific A_j(μ), m)
        # p_site now is a SINGLE VALUE, not a vector over all A
        A_j = A_known[site]             # 0..15 index
        p_site = get_p_site_for_A(site, m, bc, mat, A_j)

        p_site = max(p_site, err_correc)
        log_lik += np.log(p_site)

    return log_lik

def logliki_self_pair(parm, bc, A_known):
    """
    Log-likelihood when haplotype pair μ is known.
    parm = [logit1(epsilon)]
    bc = base counts matrix
    A_known = array of known A values for each site (0..15 index)
    A_known is sequences_pairs_mat for the specific haplotype pair being considered
    """
    m = 1
    epsilon = get_inv_logit(logit_epsilon)
    
    mat = get_p_b_diploid(m, epsilon)

    log_lik = 0.0

    for site in range(bc.shape[0]):
        # compute P(B_site | this specific A_j(μ), m)
        # p_site now is a SINGLE VALUE, not a vector over all A
        A_j = A_known[site]             # 0..15 index
        p_site = get_p_site_for_A(site, m, bc, mat, A_j)

        p_site = max(p_site, err_correc)
        log_lik += np.log(p_site)

    return log_lik

if __name__ == "__main__":
    
    ### SETUP HAPLOTYPE MATCHING DATABASE ###

    # OPTIMIZED VERSION
    
    # read sequences 
    seq_lines = [ln.strip().split() for ln in open(sequences_file) if ln.strip()]
    seq_dict = {s[0]: s[1] for s in seq_lines}

    ### ANALYSE THE DATA ###
    # set up output header
    with open(outlog, "w") as f:
        f.write("id\thaplotype_pair:p_pair:ML_mu:ML_epsilon\tcontaminated\n")
    
    #  DIPLOID lookup for O(1) mapping  
    DIPLOID_lookup = {val: idx for idx, val in enumerate(DIPLOID)}
    pairs_dict = load_pickle_dict(pairs_dict_file, compress)
    n_samples = len(pairs_dict.keys())
    variable_sites_per_sample_dict = load_pickle_dict(variable_sites_per_sample_dict_path, compress)
    normalized_self_pair_val_list = []
    full_epsilon_list = []
    print(f"n_samples in pairs_dict: {n_samples}", flush=True)

    # id, estimate of m, d, deviance statistic, lower bound m, upper bound, lower bound on d, upper bound, haplotype pairs if mixed
    # loop over each file containing base count data
    for id_ in tqdm(files, desc="Processing files", unit="file", mininterval=10):
        sample_id = Path(id_).stem.replace("_base_counts", "")
        sample_variable_sites = variable_sites_per_sample_dict.get(sample_id, set())
        sample_seq = seq_dict.get(sample_id, "")
        reconstructed_sample_seq = np.array([sample_seq[i-1] for i in sample_variable_sites]) # variable sites are 1-based

        # set up file path
        fpath = os.path.join(path, id_)
        # check file is present
        if not os.path.exists(fpath):
            print(f"Missing file for id {sample_id}, skipping....", flush=True)
            continue
        # read in base count data
        base_count = np.loadtxt(fpath, delimiter='\t', skiprows=1, dtype=int)
        # if any of the sites have a base count of zero skip the file
        if (base_count[:, 1:5].sum(axis=1).min() == 0):
            print(f"Base counts of 0 for at least one site in id {sample_id}, skipping....", flush=True)
            continue
        # ensure base counts are sorted in site order in genome
        base_count = base_count[base_count[:, 0].argsort()]
        # obtain initial estimates of m
        sites = len(base_count)
        # estimate for starting value of m from mean across heterozygous sites
        m_init = 0.999
        het_mask = (base_count[:, 1:5] == 0).sum(axis=1) < 3
        # het_mask is a boolean array indicating which sites are heterozygous
        if het_mask.sum() > 0:
            bc_het = base_count[het_mask][:, 1:5]
            # bc_het is now a numpy array of base counts at heterozygous sites
            if bc_het.shape[0] == 1:
                # only one heterozygous site
                m_init = bc_het.max() / bc_het.sum()
                # m_init is the proportion of the major base at that site (cons_depth / clean_depth in Viridian)
            else:
                bc_het_max = bc_het.max(axis=1)
                bc_het_tot = bc_het.sum(axis=1)
                m_init = np.mean(bc_het_max / bc_het_tot)
                # m_init is the mean proportion of the major base across heterozygous sites
                # (average proportion of major sample or 1 - average proportion of minor sample (which can be the contaminant))
        # avoid exact 50/50 mix
        if m_init == 0.5:
            m_init = 0.50001

        if m_init < 0.5:
            print("This sample has mu < 0.5 (contaminant supposedly in the majority)", sample_id)
            break
        #base error probability, probability that any one base call represents an error
        epsilon_init = 2e-3 # 0.2%

        parm = [get_logit2(m_init), get_logit(epsilon_init)]
        # ranges between ]0.5 ; 1.0[ and [0.04% ; 10%]
        bounds = [(get_logit2(0.500001), get_logit2(.99999)), (get_logit(1e-4), get_logit(10e-2))]
        sample_seqs = []
        sequences_pairs = []
        contaminant_list = []
        sequences_pairs_liki = []
        md_ML_m_list = []
        md_ML_epsilon_list = []
        sample_cont_list = pairs_dict.get(sample_id, [])
        n_cont = len(sample_cont_list) - 1 # Don't count self pair
        pi_c = 1/n_cont
        for contaminant in sample_cont_list:
            pair = f"{sample_id}x{contaminant}"
            cont_seq = seq_dict.get(contaminant, "")
            # Only keep variable sites for this sample
            reconstructed_cont_seq = np.array([cont_seq[i-1] for i in sample_variable_sites]) # variable sites are 1-based
            sample_seqs.append(reconstructed_cont_seq)
            sequences_pairs.append(f"{sample_id}x{contaminant}")
            contaminant_list.append(contaminant)
            # Build DIPLOID sequences_pairs_mat for this sample
            A_known = []
            for i, base in enumerate(reconstructed_sample_seq):
                for j, cont_base in enumerate(reconstructed_cont_seq):
                    diploid = f"{base}{cont_base}"
                    diploid_idx = DIPLOID_lookup[diploid]
                    A_known.append(diploid_idx)
            A_known = np.array(A_known)

            # Compute mu for the pair
            if contaminant == sample_id:
                parm = [parm[1]]
                bounds = [bounds[1]]
                def neg_logliki_epsilon(parm):
                    return -logliki_self_pair(parm, bc=base_count, A_known=A_known)
                
                opt_mu = minimize(neg_logliki_epsilon,
                                    parm,
                                    method="L-BFGS-B",
                                    bounds=bounds)
                md_ML_m = 1
                md_ML_epsilon = get_inv_logit(opt_mu.x[0])
            
            else:
                def neg_logliki_mu_epsilon(parm):
                    return -logliki_m_epsilon(parm, bc=base_count, A_known=A_known)

                opt_mu = minimize(neg_logliki_mu_epsilon,
                                parm,
                                method="L-BFGS-B",
                                bounds=bounds)
                
                md_ML_m = get_inv_logit2(opt_mu.x[0])
                md_ML_epsilon = get_inv_logit(opt_mu.x[1])

            mat_ML = get_p_b_diploid(md_ML_m, md_ML_epsilon)
            bc_liki = np.array([get_p_site(i, md_ML_m, base_count, mat_ML)
                                for i in range(base_count.shape[0])])
            pair_liki = np.sum(np.log([
                bc_liki[i, int(A_known[i])] + err_correc
                for i in range(bc_liki.shape[0]) # over all sites
            ]))
            
            if contaminant == sample_id:
                self_pair_liki = pair_liki
                self_md_ML_m = md_ML_m
                self_md_ML_epsilon = md_ML_epsilon

            else:
                pair_liki = pi_c * pair_liki
                sequences_pairs_liki.append(pair_liki)
                md_ML_m_list.append(md_ML_m)
                md_ML_epsilon_list.append(md_ML_epsilon)

        sequences_pairs_pre_map = np.array(sequences_pairs_liki)
        # We assume uniform prior over the contaminants (except self_pair)
        sequences_pairs_pre_map = sequences_pairs_pre_map + np.log(pi_con)
        # sequences_pairs_pre_map is now log posterior (up to a constant) for each pair
        mx = -max(sequences_pairs_pre_map.max(), self_pair_liki)
        # mx is used to avoid numerical underflow
        sequences_pairs_map = np.exp(sequences_pairs_pre_map + mx)
        self_pair_map = np.exp(self_pair_liki + mx)
        # sequences_pairs_map is now unnormalized posterior for each pair
        normalized_sequences_pairs_map = sequences_pairs_map / (sequences_pairs_map.sum() + self_pair_map)
        normalized_self_pair_map = self_pair_map / (sequences_pairs_map.sum() + self_pair_map)
        # normalize to get posterior probabilities
        sequences_pairs_output = f"self_pair:{normalized_self_pair_map:.3f}:{self_md_ML_m:.3f}:{self_md_ML_epsilon:.3f}"

        for idx in range(len(sequences_pairs_liki)):
            val = normalized_sequences_pairs_map[idx]
            # posterior probability for this pair
            mu_pair = md_ML_m_list[idx]
            # likelihood opt mu for this pair
            epsilon_pair = md_ML_epsilon_list[idx]
            # likelihood opt epsilon for this pair
            name = sequences_pairs[idx]
            sequences_pairs_output += f", {name}:{val:.3f}:{mu_pair:.3f}:{epsilon_pair:.3f}"

        P_cont = (pi_con * np.sum(sequences_pairs_map)) / ((1 - pi_con)*self_pair_map + (pi_con * np.sum(sequences_pairs_map)))
    
        with open(outlog, "a") as f:
            f.write(
                f"{sample_id}\t{sequences_pairs_output.strip(', ')}\t{P_cont}\n"
            )

    # print a hist plot of normalized_self_pair_val_list
    # clip normalized self pair val list between 0 and 1
    # normalized_self_pair_val_list = np.clip(normalized_self_pair_val_list, 0, 1)
    # val_hist = ff.create_distplot(
    #     [normalized_self_pair_val_list],
    #     ['Normalized self pair value'],
    #     histnorm='probability density',
    #     show_curve=True,
    #     show_hist=True,
    #     bin_size=0.05,
    #     show_rug=True,
    # )
    
    # val_hist.write_html(f"{data_dir}/8/eyre_model_{masked_or_random}_self_pair_values.html")
    # val_hist.write_image(f"{data_dir}/8/eyre_model_{masked_or_random}_self_pair_values.png")
    
    # epsilon_hist = ff.create_distplot(
    #     [full_epsilon_list],
    #     ['Error probability'],
    #     histnorm='probability density',
    #     show_curve=True,
    #     show_hist=True,
    #     bin_size=0.005,
    #     show_rug=True,
    # )
    
    # epsilon_hist.write_html(f"{data_dir}/8/eyre_model_{masked_or_random}_epsilon.html")
    # epsilon_hist.write_image(f"{data_dir}/8/eyre_model_{masked_or_random}_epsilon.png")