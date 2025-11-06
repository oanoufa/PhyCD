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
from scipy.special import expit, logit
from itertools import combinations, product
from collections import defaultdict
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd
import random
import os
from tqdm import tqdm
from pathlib import Path
from time import time

### OPTIONS ###

#path to directory with base counts
path = Path("/nfs/research/goldman/anoufa/src/test_eyre_model/pango_cons_seq_lev_4/")
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

#file to record output
outlog = path / "sample_output.txt"

do_bootstrap = False  #whether to perform bootstrap to obtain confidence intervals
t0 = time()

# epsilon to use in place of 0 probabilities to avoid log(0)
eps = 1e-300

#number of cores to use for bootstrapping on Mac / Linux, ignored on Windows
bs_cores = 4
#number of bootstrap iterations to run for parameter estimates
bs_iter = 1000

#base error probability, probability that any one base call represents an error
p_err = 2e-3

#deviance statistic threshold, above this perform haplotype matching, below assume not mixed
dev_thres = 19.4

### END OPTIONS ###

### CONSTANTS ###
DIPLOID = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
BASES = ["A", "C", "G", "T"]


### FUNCTIONS ###
#functions for defining m, mixture proportion and d, between haplotype diversity by maximum likelihood

#define functions to keep values between 0.5 and 1
def get_logit2(p):
    return math.log((2*p-1)/(1-(2*p-1)))

def get_inv_logit2(p):
    return 0.5 + (math.exp(p) / (1 + math.exp(p))) / 2


#define functions to keep values between 0 and 1
def get_logit(p):
    return math.log(p/(1-p))
def get_inv_logit(p):
    return math.exp(p) / (1 + math.exp(p))

#several functions to generate the likelihood of m, the mixture proportion
def get_p_b(A, m, b, e):
    # returns p(b|A,m,e) where A = {AA, AC, ..., TT} and b is a base in a single read, and e is the read error probability
    p_b_1 = m * e/3  # p where major ≠ b 
    p_b_2 = (1-m) * e/3  # p where minor ≠ b
    if b == A[0]:
        p_b_1 = m * (1-e)  # p where major = b
    if b == A[1]:
        p_b_2 = (1-m) * (1-e)  # p where minor = b
    p_b = p_b_1 + p_b_2
    return p_b

def get_p_b_diploid(m):
    # returns a matrix of p(b|A,m,e) for all values of b {A, C, G, T} - rows and values of A - columns
    # e is set by constant above
    mat = np.empty((16, 4))
    for i, A in enumerate(DIPLOID):
        mat[i, :] = [get_p_b(A, m, b, p_err) for b in BASES]
    return mat

def get_p_site(site, m, bc, mat):
    # returns p(B1 to Bn | A, m, e) where B1 to Bn are all reads at a given site
    # - product of p(b|A,m,e) across all reads - taken from base counts
    # where mat is a 16x4 matrix from get.p.b.diploid
    counts = bc[site, 1:5]  # ensure four BASES
    return np.exp(np.sum(np.log(mat) * counts, axis=1))
    
def get_d(d):
    # returns a vector of probabilities for the 16 DIPLOID forms
    aa = (1-d)/4
    ab = d/12
    return np.array([aa, ab, ab, ab, ab, aa, ab, ab, ab, ab, aa, ab, ab, ab, ab, aa])

def logliki_m(parm, bc):
    # return the log likihood of m,d - product over all sites [ sum over all values of A [ p(B1-Bn | A, m, e)  * p(A) ] ]
    # p(A) from the vector given by get.d
    m = get_inv_logit2(parm[0])
    d = get_inv_logit(parm[1])
    mat = get_p_b_diploid(m)
    probs_A = get_d(d)
    log_lik = 0.0
    eps = 1e-300  # to avoid log(0)

    for site in range(bc.shape[0]):
        p_site = get_p_site(site, m, bc, mat)
        p_site = np.maximum(p_site, eps)
        # use log-sum-exp trick for numerical stability
        log_lik += np.logaddexp.reduce(np.log(p_site) + np.log(probs_A))
    return log_lik

def get_dev(ml, bc):
    #compare the likelihood obtained to the likelihood under an approximation to the null hypothesis
    h0_liki = logliki_m([get_logit2(0.9999), get_logit(0.00001)], bc)
    D = -2 * (h0_liki - ml)
    return D

def sample_bc(bc):
    #function for sampling base counts
    np.random.seed(42)
    idx = np.random.choice(bc.shape[0], bc.shape[0], replace=False)
    return bc[idx]


if __name__ == "__main__":
    
    ### SETUP HAPLOTYPE MATCHING DATABASE ###
    print(f"Loading haplotype database... t0: {(time() - t0):.2f}", flush=True)
    # # read in haplotypes
    # sequences = pd.read_table(sequences_file, header=None)

    # # read in population frequency of each haplotype
    # sequences_freq = pd.read_table(sequences_freq_file, header=None, index_col=0)

    # # find the least common haplotype
    # min_freq = sequences_freq.iloc[:, 0].min()

    # # split sequence strings into matrix of sites
    # sequences_mat = np.array([list(seq) for seq in sequences.iloc[:, 1]])
    # sequences_mat = pd.DataFrame(sequences_mat, index=sequences.iloc[:, 0])

    # sequences_n = sequences_mat.shape[0]   # number of haplotypes
    # sequences_len = sequences_mat.shape[1] # number of variable sites
    # sequences_comb = sequences_n * (sequences_n - 1)  # number of combinations of different minor / major haplotypes

    # # generate vector to contain all haplotype pairs
    # sequences_pairs = []
    # # generate vector of relative frequencies of different haplotype pairs
    # sequences_pairs_freq = []

    # # populate vector of haplotype pair names and generate vector of relative frequencies of different haplotype pairs
    # for i, j in product(sequences_mat.index, repeat=2):
    #     if i != j:
    #         sequences_pairs.append(f"{i}x{j}")
    #         freq_i = sequences_freq.loc[i, sequences_freq.columns[0]] if i in sequences_freq.index else min_freq
    #         freq_j = sequences_freq.loc[j, sequences_freq.columns[0]] if j in sequences_freq.index else min_freq
    #         sequences_pairs_freq.append(freq_i * freq_j)

    # # normalise pair frequency data to give prior probability of each pair
    # sequences_pairs_freq = np.array(sequences_pairs_freq)
    # sequences_pairs_prior = sequences_pairs_freq / sequences_pairs_freq.sum()

    # # generate matrix of sequences for all haplotype pairs
    # sequences_pairs_mat = pd.DataFrame(
    #     index=sequences_pairs,
    #     columns=range(sequences_len),
    #     dtype=float
    # )

    # for i, j in product(sequences_mat.index, repeat=2):
    #     if i != j:
    #         sequences_pairs_mat.loc[f"{i}x{j}"] = [
    #             DIPLOID.index(f"{a}{b}") + 1  # +1 to match R's 1-based index
    #             for a, b in zip(sequences_mat.loc[i], sequences_mat.loc[j])
    #         ]
    
    
    # OPTIMIZED VERSION
    
    # --- read sequences (fast, no pandas) ---
    seq_lines = [ln.strip().split() for ln in open(sequences_file) if ln.strip()]
    seq_ids = [s[0] for s in seq_lines]
    seq_strs = [s[1] for s in seq_lines]

    print(f"Read {len(seq_ids)} haplotypes. t1: {(time() - t0):.2f}", flush=True)

    # --- read frequencies into dict (fallback to minimal later) ---
    freq_lines = [ln.strip().split() for ln in open(sequences_freq_file) if ln.strip()]
    freq_dict = {k: float(v) for k, v in freq_lines}
    min_freq = min(freq_dict.values()) if freq_dict else 1.0

    print(f"Read frequencies for {len(freq_dict)} haplotypes. t2: {(time() - t0):.2f}", flush=True)

    # --- convert to numpy char matrix ---
    sequences_mat = np.array([list(seq) for seq in seq_strs], dtype="<U1")  # shape (n, L)
    sequences_n, sequences_len = sequences_mat.shape
    
    print(f"Converted sequences to matrix. t3: {(time() - t0):.2f}", flush=True)

    # --- number of ordered pairs (i != j) and generate pairs list ---
    sequences_comb = sequences_n * (sequences_n - 1)
    pairs = [(seq_ids[i], seq_ids[j]) for i in range(sequences_n) for j in range(sequences_n) if i != j]
    sequences_pairs = [f"{a}x{b}" for a, b in pairs]
    
    print(f"Generated {sequences_comb} haplotype pairs. t4: {(time() - t0):.2f}", flush=True)

    # --- compute pair priors (vectorized comprehension) ---
    pair_freqs = np.array([freq_dict.get(a, min_freq) * freq_dict.get(b, min_freq) for a, b in pairs], dtype=float)
    sequences_pairs_prior = pair_freqs / pair_freqs.sum()

    print(f"Computed pair priors. t5: {(time() - t0):.2f}", flush=True)

    # --- DIPLOID lookup for O(1) mapping  ---
    DIPLOID_lookup = {val: idx for idx, val in enumerate(DIPLOID)}

    print(f"Prepared DIPLOID lookup. t6: {(time() - t0):.2f}", flush=True)

    # --- build sequences_pairs_mat as a numpy array (rows = pair, cols = sites) ---
    sequences_pairs_mat = np.empty((sequences_comb, sequences_len), dtype=np.int16)
    # Precompute index mapping from seq_id to row index in sequences_mat
    id_to_idx = {seq_ids[i]: i for i in range(sequences_n)}

    # Precompute lookup table as 4x4 array (if DIPLOID_lookup has simple mapping)
    lut = np.zeros((4, 4), dtype=np.uint8)
    for i, x in enumerate(BASES):
        for j, y in enumerate(BASES):
            lut[i, j] = DIPLOID_lookup[f"{x}{y}"]

    # Map sequences_mat characters to indices 0-3
    base_to_idx = {b: i for i, b in enumerate(BASES)}
    seq_idx = np.vectorize(base_to_idx.get)(sequences_mat)

    # For each pair, just lookup
    for row_idx, (a, b) in enumerate(tqdm(pairs, desc="Building sequences_pairs_mat", unit="pair", mininterval=10)):
        ia, ib = id_to_idx[a], id_to_idx[b]
        sequences_pairs_mat[row_idx, :] = lut[seq_idx[ia], seq_idx[ib]]


    print(f"Haplotype database loaded. t7: {(time() - t0):.2f}", flush=True)

    ### ANALYSE THE DATA ###
    # set up output header
    with open(outlog, "w") as f:
        if do_bootstrap:
            f.write("id\tML_m\tML_d\tdeviance\tm0.025\tm0.975\td0.025\td0.975\tp_haplotype_pair\n")
        else:
            f.write("id\tML_m\tML_d\tdeviance\tp_haplotype_pair\n")
    # id, estimate of m, d, deviance statistic, lower bound m, upper bound, lower bound on d, upper bound, haplotype pairs if mixed

    # loop over each file containing base count data
    for id_ in tqdm(files, desc="Processing files", unit="file"):
        # set up file path
        fpath = os.path.join(path, id_)
        # check file is present
        if not os.path.exists(fpath):
            print(f"Missing file for id {id_}, skipping....", flush=True)
            continue

        print(f"Starting {id_}.... time: {(time() - t0):.2f}", flush=True)

        # read in base count data
        base_count = pd.read_csv(fpath, sep="\t")

        # if any of the sites have a base count of zero skip the file
        if (base_count.iloc[:, 1:5].sum(axis=1).min() == 0):
            print(f"Base counts of 0 for at least one site in id {id_}, skipping....", flush=True)
            continue

        # ensure base counts are sorted in site order in genome
        base_count = base_count.sort_values(by=base_count.columns[0]).reset_index(drop=True)

        # obtain initial estimates of d and m
        sites = len(base_count)
        d_obs = np.sum((base_count.iloc[:, 1:5] == 0).sum(axis=1) < 3)
        if d_obs == 0:
            d_obs = 1e-5
        d_init = d_obs / sites

        # estimate for starting value of m from mean across heterozygous sites
        m_init = 0.999
        het_mask = (base_count.iloc[:, 1:5] == 0).sum(axis=1) < 3
        if het_mask.sum() > 0:
            bc_het = base_count.loc[het_mask, base_count.columns[1:5]].to_numpy()
            if bc_het.shape[0] == 1:
                m_init = bc_het.max() / bc_het.sum()
            else:
                bc_het_max = bc_het.max(axis=1)
                bc_het_tot = bc_het.sum(axis=1)
                m_init = np.mean(bc_het_max / bc_het_tot)

        # avoid exact 50/50 mix
        if m_init == 0.5:
            m_init = 0.50001

        init_parm = [get_logit2(m_init), get_logit(d_init)]

        ### numerically optimise the values of m and d and obtain ML
        def neg_logliki(par):
            return -logliki_m(par, base_count.to_numpy())

        opt = minimize(neg_logliki, init_parm, method="BFGS")
        md_ML_m = get_inv_logit2(opt.x[0])
        md_ML_d = get_inv_logit(opt.x[1])

        md_dev = get_dev(opt.fun * -1, base_count.to_numpy())

        print(f"got ML m: {md_ML_m} and d: {md_ML_d} (deviance: {md_dev}). time: {(time() - t0):.2f}", flush=True)

        # run bootstrap
        def bootstrap_iter(_):
            bc_samp = sample_bc(base_count.to_numpy())
            opt_bs = minimize(lambda p: -logliki_m(p, bc_samp), init_parm, method="BFGS")
            bs_m = get_inv_logit2(opt_bs.x[0])
            bs_d = get_inv_logit(opt_bs.x[1])
            return [bs_m, bs_d]

        if do_bootstrap:
            if os.name != "nt":
                bs_md = Parallel(n_jobs=bs_cores)(delayed(bootstrap_iter)(i) for i in range(bs_iter))
            else:
                bs_md = [bootstrap_iter(i) for i in range(bs_iter)]
            bs_md = np.array(bs_md)

            bs_m_ci = np.quantile(bs_md[:, 0], [0.025, 0.975])
            bs_d_ci = np.quantile(bs_md[:, 1], [0.025, 0.975])

            print(f"m ({bs_m_ci[0]} - {bs_m_ci[1]})")
            print(f"d ({bs_d_ci[0]} - {bs_d_ci[1]})")

        # if reaches deviance threshold then perform haplotype matching
        sequences_pairs_output = ""
        if md_dev > dev_thres:
            print(f"Performing haplotype matching for {id_}.... time: {(time() - t0):.2f}", flush=True)
            
            mat_ML = get_p_b_diploid(md_ML_m)
            bc_liki = np.array([get_p_site(i, md_ML_m, base_count.to_numpy(), mat_ML)
                                for i in range(base_count.shape[0])])

            sequences_pairs_liki = np.array([
                np.sum(np.log([
                    bc_liki[i, int(sequences_pairs_mat[x, i])] + eps
                    for i in range(bc_liki.shape[0])
                ]))
                for x in range(sequences_comb)
            ])

            print(f"Computed sequence pair likelihoods. time: {(time() - t0):.2f}", flush=True)

            sequences_pairs_pre_map = sequences_pairs_liki + np.log(sequences_pairs_prior)
            mx = -np.max(sequences_pairs_pre_map)
            sequences_pairs_map = np.exp(sequences_pairs_pre_map + mx)
            sequences_pairs_map /= sequences_pairs_map.sum()

            top_idx = np.argmax(sequences_pairs_map)
            print(sequences_pairs[top_idx], flush=True)

            sorted_idx = np.argsort(-sequences_pairs_map)
            cum_sum = np.cumsum(sequences_pairs_map[sorted_idx])

            for k, idx in enumerate(sorted_idx):
                if k == 0 or cum_sum[k - 1] <= 0.99:
                    val = sequences_pairs_map[idx]
                    name = sequences_pairs[idx]
                    print(val)
                    sequences_pairs_output += f", {name}:{val:.5f}"

        with open(outlog, "a") as f:
            if do_bootstrap:
                f.write(
                    f"{id_}\t{md_ML_m}\t{md_ML_d}\t{md_dev}\t"
                    f"{bs_m_ci[0]}\t{bs_m_ci[1]}\t{bs_d_ci[0]}\t{bs_d_ci[1]}\t"
                    f"{sequences_pairs_output.strip(', ')}\n"
                )
            else:
                f.write(
                    f"{id_}\t{md_ML_m}\t{md_ML_d}\t{md_dev}\t"
                    f"{sequences_pairs_output.strip(', ')}\n"
                )