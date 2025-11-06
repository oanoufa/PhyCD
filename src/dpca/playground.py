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

def build_maple_vcf_dict(clean_file_path, param_term, masked_or_random):
    """Build a dictionary mapping sample names and type (random, masked, unmasked) to their MAPLE VCF entries."""
    maple_vcf_dict = {}
    # FIRST RETRIEVE ALL SAMPLES FROM CLEAN TREE
    print("Building MAPLE VCF dictionary from clean tree alignment...")
    with open(clean_file_path, "rt") as file:
        file.readline()  # SKIP REFERENCE
        line = file.readline().strip()
        while line:
            if line.startswith(">"):
                sample_name = line[1:]  # Remove '>'
                maple_vcf_dict[(sample_name, "clean")] = []
                line = file.readline().strip()
                while line and not line.startswith(">"):
                    maple_vcf_dict[(sample_name, "clean")].append(line)
                    line = file.readline().strip()
            else:
                line = file.readline().strip()

    alignment_folder = Path("/nfs/research/goldman/anoufa/data/MAPLE_input/")
    print(f"Adding unmasked and {masked_or_random} entries...")
    for alignment_file in tqdm(
        alignment_folder.glob(f"*_alignment_{param_term}.maple"),
        desc="Building MAPLE VCF dictionary",
        mininterval=60,
    ):
        # OPEN THE 3 FILES unmasked_alignment, masked_alignment, random_alignment
        # Typical file name /nfs/research/goldman/anoufa/data/MAPLE_input/DPCA_unmasked_alignment_1_0.1_1_0_30000.maple
        stem = alignment_file.stem
        if "unmasked" in stem:
            file_type = "unmasked"
        elif "masked" in stem:
            file_type = "masked"
            if masked_or_random != "masked":
                continue
        elif "random" in stem:
            file_type = "random"
            if masked_or_random != "random":
                continue
        else:
            raise ValueError(f"Unknown file type in {alignment_file.name}")

        with open(alignment_file, "rt") as file:
            # Skip reference (first entry)
            file.readline()
            line = file.readline().strip()

            while line:
                if line.startswith(">"):
                    entry_name = line[1:]  # Remove '>'
                    sample_name = entry_name.split("_", 1)[
                        1
                    ]  # Get the sample name after the first underscore
                    maple_vcf_dict[(sample_name, file_type)] = []
                    line = file.readline().strip()
                    while line and not line.startswith(">"):
                        maple_vcf_dict[(sample_name, file_type)].append(line)
                        line = file.readline().strip()
                else:
                    line = file.readline().strip()

    return maple_vcf_dict


def look_for_maple_vcf(sample_to_look_for, cons_or_masked_or_random, clean=False):
    
    # Define the path to the folder containing the tsv files
    folder_path = Path("/nfs/research/goldman/anoufa/data/dpca/batches")    
    maple_entry = []

    # Find the sample's position in the batch file
    path_vdn = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"
    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        for i, line in enumerate(f):
            if i == 0:
                continue
            # Get the read_name and path
            read_name, path = line.strip().split("\t")
            if read_name == sample_to_look_for:
                # Retrieve batch number
                batch_to_look_for = i % n_batch

    if clean:
        file_to_look_for = f"/nfs/research/goldman/anoufa/data/MAPLE_input/clean_tree_alignment_file_{param_term}.maple"
        # Read the tsv file
        with open(file_to_look_for, "rt") as file:
        
            # Extract the sample lengths for the three types of sequences
            for line in file:
                # Skip first line (header) and empty lines
                line = line.strip()
                entry_name = f">{sample_to_look_for}"
                if line == entry_name:
                    line = next(file).strip()
                    while not line.startswith(">"):
                        maple_entry.append(line)
                        line = next(file).strip()
                    return maple_entry
                
    else:
        file_to_look_for = f"{folder_path}/maple_alignment_batch{batch_to_look_for}_{param_term}.gz"

        # Read the tsv file
        with gzip.open(file_to_look_for, "rt") as file:
        
            # Extract the sample lengths for the three types of sequences
            for line in file:
                # Skip first line (header) and empty lines
                line = line.strip()

                if line == f">{cons_or_masked_or_random}_{sample_to_look_for}":
                    line = next(file).strip()
                    while not line.startswith(">"):
                        maple_entry.append(line)
                        line = next(file).strip()
                    return maple_entry
                
    print(f"Could not find {cons_or_masked_or_random} entry for {sample_to_look_for} in {file_to_look_for}")
    return False


def load_globals():
    """Load global dictionaries from pickle files or build them if they don't exist.
    """
    
    global mutation_to_variants, maple_vcf_dict

    # Load the mutation to variants dictionary once
    mutation_to_variants_path = Path(
        f"/nfs/research/goldman/anoufa/data/MAPLE_input/dict/clean_samples_mut_to_var_dict_{param_term}.pickle"
    )

    clean_tree_path_1 = Path(f"/nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_tree_alignment_file_{param_term}.maple")

    if clean_tree_path_1.exists():
        clean_tree_alignment_path = clean_tree_path_1
    else:
        clean_tree_alignment_path = Path(
            f"/nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/clean_tree_alignment_file_{param_term}.maple"
        )

    
    maple_vcf_dict_path = Path(
        f"/nfs/research/goldman/anoufa/data/MAPLE_input/dict/maple_vcf_dict_{param_term}_{masked_or_random}.pickle"
    )
    
    # IF THE FILE EXISTS LOAD IT, OTHERWISE BUILD IT AND STORE IT
    if maple_vcf_dict_path.exists():
        print(f"Loading existing MAPLE VCF dictionary from {maple_vcf_dict_path}...", flush=True)
        with open(maple_vcf_dict_path, "rb") as f:
            maple_vcf_dict = pickle.load(f)
        print(f"Loaded existing MAPLE VCF dictionary with {len(maple_vcf_dict)} entries.", flush=True)
    else:
        maple_vcf_dict = build_maple_vcf_dict(
            clean_tree_alignment_path, param_term, masked_or_random
        )
        print("Storing MAPLE VCF dictionary...", flush=True)
        with open(maple_vcf_dict_path, "wb") as f:
            pickle.dump(maple_vcf_dict, f)
        print(f"Stored MAPLE VCF dictionary with {len(maple_vcf_dict)} entries to {maple_vcf_dict_path}.", flush=True)
    
    if mutation_to_variants_path.exists():
        with open(mutation_to_variants_path, "rb") as f:
            mutation_to_variants = pickle.load(f)
        print(f"Loaded existing mutation to variants dictionary with {len(mutation_to_variants)} mutations.", flush=True)
    else:
        print('Could not find mutation to variants dictionary, please rerun "2_process_gmf_output.py".', flush=True)

def look_for_closest_variant_optimized_unified(
    input_tsv,
    output_tsv,
    mutation_col_name,
    masked_or_random,
    ref_seq,
):
    """
    Iterate over the input TSV file, for each sample find the closest variants in the clean tree
    based on the mutations present in the sample and the mutations associated with each variant.
    On those closest variants, compute the hamming distance to the unmasked sequence at the masked regions.
    Write the results to the output TSV file with additional columns:
    - closest_variant: the name(s) of the closest variant(s) in the clean tree
    - n_candidates: number of candidate variants with the highest closeness score
    - closeness: the highest closeness score (number of shared mutations)
    - closeness_ratio: closeness / number of mutations in the sample
    - hamming_dist: the best hamming distance (proportion of matches) among the closest variants
    """

    print(f"[PID {os.getpid()}] Worker started with {input_tsv}", flush=True)
    # LOAD THE DICTIONARIES
    print(f"[PID {os.getpid()}] Loaded mutation to variants dictionary with {len(mutation_to_variants)} mutations.", flush=True)
    print(f"[PID {os.getpid()}] Loaded MAPLE VCF dictionary with {len(maple_vcf_dict)} entries.", flush=True)
    
    # Main processing
    with open(input_tsv, "r") as f_in, open(output_tsv, "w", newline="") as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")

        # Write new header
        header = next(reader)
        sample_idx = header.index("sample_name")
        mutation_idx = header.index(mutation_col_name)

        # Prepare output header
        new_header = header + [
            "closest_variant",
            "n_candidates",
            "closeness",
            "closeness_ratio",
        ]

        writer.writerow(new_header)

        for cols in tqdm(reader, desc=f"[PID {os.getpid()}] Processing samples", mininterval=120):
            sample_name = cols[sample_idx]
            sample_mut = cols[mutation_idx]

            # Default output values
            closest_variant = "No mutations masked by gen_maple_file"
            n_candidates = 0
            closeness = 0
            closeness_ratio = 0
            
            # Retrieve VCF entry of the sample, only keep the entries with SNPs, add entries for heterozygous sites
            unmasked_vcf_entry = maple_vcf_dict.get(
                (sample_name, "unmasked"), []
            )
            
            masked_vcf_entry = maple_vcf_dict.get(
                (sample_name, masked_or_random), []
            )
            
            # Parse both entries to only keep the regions masked by our process
            masked_regions_masked = parse_intervals(
                masked_vcf_entry
            )  # MASKED REGIONS IN THE MASKED OR RANDOM ENTRY
            masked_regions_unmasked = parse_intervals(
                unmasked_vcf_entry
            )  # MASKED REGIONS IN THE UNMASKED (VIRIDIAN) ENTRY
            
            masked_regions_alg = subtract_intervals(
                masked_regions_masked, masked_regions_unmasked
            )  # REGIONS MASKED BY OUR PROCESS ONLY
            
            # Only keep the VCF entries of unmasked_vcf_entry that fall within masked_regions_alg
            alg_masked_and_het_mut = set()
            
            for mut in parse_mutations(unmasked_vcf_entry):
                if len(mut) == 3:
                    # Entry is a deletion or Ns, we don't consider it. example: (n, 8765, 14)
                    # Maybe we should consider deletions?
                    continue
                else:
                    # Entry is a SNP, example: (t, 2343)
                    # Format it to a string like 2343T
                    base, pos = mut
                    str_mut = f"{pos}{base}".upper()
                    for start, end in masked_regions_alg:
                        if start <= pos <= end:
                            alg_masked_and_het_mut.add(str_mut)
                            break
            # Add heterozygous sites in the entries
            sample_het_sites = het_sites_dict.get(sample_name, [])
            for pos, minor_allele, _ in sample_het_sites:
                str_het_mut = f"{pos}{minor_allele}".upper()
                alg_masked_and_het_mut.add(str_het_mut)

            # alg_masked_vcf_entry now contains only the mutations in the unmasked entry that fall within the regions masked by our process, plus the heterozygous sites
            # This is everything we want to consider for finding the closest variants

            # Expand mutations masked by gen_maple_file
            expanded_mut_set = set()
            for mut in alg_masked_and_het_mut:
                expanded_mut_set.update(expand_ambiguous_mutation(mut, remove_starting_nt=False))

            # Find candidate variants
            candidate_counts = Counter()
            for mut in expanded_mut_set:
                candidate_counts.update(mutation_to_variants.get(mut, []))

            if candidate_counts:
                max_closeness = max(candidate_counts.values())
                best_contaminant = [
                    v
                    for v, count in candidate_counts.items()
                    if count == max_closeness
                ]

                # Save computed values
                closest_variant = best_contaminant
                n_candidates = len(best_contaminant)
                closeness = max_closeness
                closeness_ratio = max_closeness / len(sample_mut)

            # Write updated row
            writer.writerow(
                cols
                + [
                    closest_variant,
                    n_candidates,
                    closeness,
                    closeness_ratio,
                ]
            )
    print(f"Updated TSV written to: {output_tsv}", flush=True)

if __name__ == "__main__":

    # OPEN /nfs/research/goldman/anoufa/data/MAPLE_input/others/viridian_samples.metadata.tsv.gz
    # AND CHECK THE NUMBER OF SAMPLES IN IT (n_lines)
    # 5 959 032 samples (should have all my samples)
    
    # OPEN /nfs/research/goldman/anoufa/data/MAPLE_output/processed_placements/processed_placements_results_with_contaminants_random_0.1_3.4.tsv
    # CHECK WHY ITS 630Mb instead of the usual ~8Mb
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.expand_frame_repr', False)

    paths_df = [Path("/nfs/research/goldman/anoufa/data/MAPLE_output/processed_placements/processed_placements_results_with_contaminants_random_0.1_3.4.tsv"),
                Path("/nfs/research/goldman/anoufa/data/MAPLE_output/processed_placements/processed_placements_results_with_contaminants_masked_0.1_3.15.tsv")]
    
    for path_df in paths_df:
        print(f"Processing {path_df}")
        df_pp = pd.read_csv(path_df, sep="\t")

        # Print head FULL
        print(df_pp.head(5))
        print(df_pp.shape)
        
        # statistics on the df (min max mean closeness, closeness_ratio, n_candidates)
        print("Closeness statistics:")
        print(f"Min closeness: {df_pp['closeness'].min()}")
        print(f"Max closeness: {df_pp['closeness'].max()}")
        print(f"Mean closeness: {df_pp['closeness'].mean()}")
        print(f"Min closeness ratio: {df_pp['closeness_ratio'].min()}")
        print(f"Max closeness ratio: {df_pp['closeness_ratio'].max()}")
        print(f"Mean closeness ratio: {df_pp['closeness_ratio'].mean()}")
        print(f"Min n_candidates: {df_pp['n_candidates'].min()}")
        print(f"Max n_candidates: {df_pp['n_candidates'].max()}")
        print(f"Mean n_candidates: {df_pp['n_candidates'].mean()}")
        
        # Show row with max n_candidates
        max_n_candidates = df_pp['n_candidates'].max()
        print("Row(s) with max n_candidates:")
        print(len(df_pp[df_pp['n_candidates'] == max_n_candidates]['closest_variant'].values[0].split(";")))
        print("\n\n")