from pathlib import Path
import gzip
import lzma
import _params
import pickle
import csv
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
from collections import Counter
from itertools import groupby
import lzma
import gzip
import sys
sys.path.append("/nfs/research/goldman/anoufa/src/dpca")
from _aux_functions import parse_maple_file, generate_sample_list

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
path_psmp_masked = _params.processed_placement_masked
path_psmp_random = _params.processed_placement_random
param_term = _params.param_term
num_cores = _params.num_cores

if __name__ == "__main__":

    # OPEN /nfs/research/goldman/anoufa/data/MAPLE_input/others/viridian_samples.metadata.tsv.gz
    # AND CHECK THE NUMBER OF SAMPLES IN IT (n_lines)
    # 5 959 032 samples (should have all my samples)
    
    # Generate a dataset S2 to see how the Eyre model scales.
    path_to_variants_file = Path("/nfs/research/goldman/anoufa/data/MAPLE_input/others/pango-consensus-sequences_level_4_20230228.maple")
    
    # This file is MAPLE VCF file with the pango consensus sequences at level 4.
    # From this file we will extract the sites of variability (all the sites present at least in one of the sequences)
    ref_name, ref_seq, maple_dict = parse_maple_file(path_to_variants_file)
    variable_sites = set()
    for entry_lines in tqdm(maple_dict.values(), desc="Extracting variable sites from maple_file"):
        for line in entry_lines:
            parts = line.split()
            pos = int(parts[1])
            variable_sites.add(pos)
    print(f"Number of variable sites in pango-consensus-sequences_level_4_20230228: {len(variable_sites)}")
    print(f"Number of sequences in pango-consensus-sequences_level_4_20230228: {len(maple_dict)}")
    
    # Reconstruct sequences only at variable sites
    reconstructed_seqs = {}
    for seq_name, entry_lines in tqdm(maple_dict.items(), desc="Reconstructing sequences at variable sites"):
        seq = list(ref_seq)
        for line in entry_lines:
            parts = line.split()
            pos = int(parts[1])
            base = parts[0].upper()
            if base in ["A", "C", "G", "T"]:
                seq[pos - 1] = base  # MAPLE positions are 1-based
        # Keep only variable sites
        var_seq = [seq[pos - 1] for pos in sorted(variable_sites)]
        reconstructed_seqs[seq_name] = "".join(var_seq)
    
    # Write the reconstructed sequences to the haplotype_sequence.txt file, with two columns: sequence_name, sequence (unnamed)
    haplotype_sequence_path = Path("/nfs/research/goldman/anoufa/src/test_eyre_model/pango_cons_seq_lev_4/haplotype_sequence.txt")
    with haplotype_sequence_path.open("w") as f:
        for seq_name, seq in reconstructed_seqs.items():
            f.write(f"{seq_name}\t{seq}\n")
            
    # Generate the haplotype_frequency.txt file, with two columns: sequence_name, frequency
    # For simplicity, we will assign a frequency of 1/number_of_sequences to each sequence
    haplotype_frequency_path = Path("/nfs/research/goldman/anoufa/src/test_eyre_model/pango_cons_seq_lev_4/haplotype_frequency.txt")
    frequency = 1 / len(reconstructed_seqs)
    with haplotype_frequency_path.open("w") as f:
        for seq_name in reconstructed_seqs.keys():
            f.write(f"{seq_name}\t{frequency}\n")
    
    # Select some samples to test!
    
    # Clear previous samples
    output_folder = Path("/nfs/research/goldman/anoufa/src/test_eyre_model/pango_cons_seq_lev_4/")
    for file in output_folder.glob("*_base_counts.txt"):
        file.unlink()
    # Open a processed_placements_results_with_contaminants and select the 10 best candidates
    processed_placements_path = Path("/nfs/research/goldman/anoufa/data/MAPLE_output/processed_placements/processed_placements_results_with_contaminants_masked_0.1_3.1.tsv")
    df_pp = pd.read_csv(processed_placements_path, sep="\t")
    print(df_pp.head())

    # Sort by closeness and take the top 5 samples
    top_samples = df_pp.sort_values(by="closeness", ascending=False).head(5)["sample_name"].tolist()
    print(f"Top 5 samples selected for testing: {top_samples}")
    
    # Also take 5 clean samples from the clean alignment file
    path_clean_alignment = Path("/nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/clean_tree_alignment_file_0.1_3.maple")
    _, _, clean_maple_dict = parse_maple_file(path_clean_alignment)
    clean_sample_names = list(clean_maple_dict.keys())
    clean_top_samples = clean_sample_names[:5]
    print(f"5 clean samples selected for testing: {clean_top_samples}")

    top_samples.extend(clean_top_samples)
    # Go look for the qc files of those samples and write the counts of nt at variable sites only of those samples
    name_path_list = generate_sample_list(top_samples)
    
    for path, read_name in tqdm(name_path_list,
                                desc="Processing samples",
                                unit="sample",
                                mininterval=120):
        
        ref_pos_list = []
        unma_seq = []
        counts_list = []
        
        with gzip.open(path, "rt") as f:  # read as text directly
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)  # skip header
            
            for line in reader:
                if not line or len(line) < 19:
                    continue
                
                if line[9].isdigit():
                    ref_pos = int(line[0])
                    clean_depth = int(line[9])
                    
                    # minor allele calculation
                    counts = [
                        int(line[11]) + int(line[12]),  # A
                        int(line[13]) + int(line[14]),  # C
                        int(line[15]) + int(line[16]),  # G
                        int(line[17]) + int(line[18])   # T
                    ]
                    if clean_depth == 0:
                        continue
                    
                    # list like ref_pos, A_count, C_count, G_count, T_count
                    if ref_pos in variable_sites:
                        # check that at least one count > 0
                        if sum(counts) == 0:
                            # add 1 everywhere to avoid zero counts
                            counts = [1, 1, 1, 1]
                        counts_list.append([ref_pos] + counts)
                    
                else:
                    # Deletion ( - row) or rows full of dots (start/end gaps)
                    ref_pos = int(line[0])
                    cons_nt = line[4].lower()
                
                # track unmasked sequence
                if ref_pos_list and ref_pos == ref_pos_list[-1]:
                    unma_seq[-1] += cons_nt
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    
        # Now write the counts at variable sites only
        output_counts_path = Path(f"/nfs/research/goldman/anoufa/src/test_eyre_model/pango_cons_seq_lev_4/{read_name}_base_counts.txt")
        with output_counts_path.open("w") as f:
            # header
            f.write("site\tA\tC\tG\tT\n")
            for counts in counts_list:
                f.write("\t".join(map(str, counts)) + "\n")

    # Write the list of file names to filenames_list.txt
    filenames_list_path = Path("/nfs/research/goldman/anoufa/src/test_eyre_model/pango_cons_seq_lev_4/filenames_list.txt")
    with filenames_list_path.open("w") as f:
        for read_name in top_samples:
            f.write(f"{read_name}_base_counts.txt\n")