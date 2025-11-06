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
from _aux_functions import update_params_file, expand_ambiguous_mutation
from multiprocessing import Pool
from os import cpu_count
import os
from itertools import chain

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

if masked_or_random == "masked":
    het_sites_dict_path = _params.het_sites_dict_masked
else:
    het_sites_dict_path = _params.het_sites_dict_random


# Global variables
mutation_to_variants = None
maple_vcf_dict = None
het_sites_dict = None


def sorting_key(x):
    """
    Sorting key for the mutations.
    It sorts by the position of the mutation in the genome.
    Handle both cases like A4563T and A(0.000000/0.003493/0.000000/0.996507)27883T
    """
    if "(" in x:
        pos = x.split(")")[1]
        pos = pos[:-1]

    else:
        pos = x[1:-1]

    return int(pos)

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

    alignment_folder = Path("/nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/")
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


def parse_intervals(lines):
    """
    Parse a MAPLE entry to retrieve the masked regions as a list of (start, end) tuples.
    Each line is expected to be in the format: base pos [length]
    where base is 'N' or 'n' for masked regions, pos is the starting position, and length is optional (default 1).
    """
    intervals = []
    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue
        base, pos = parts[0], int(parts[1])
        if len(parts) == 3:
            length = int(parts[2])
        else:
            length = 1
        if base in ["n", "N"]:
            intervals.append((pos, pos + length - 1))
    return sorted(intervals)


def subtract_intervals(masked_intervals, unmasked_intervals):
    """
    Return parts of masked_intervals that do not overlap with unmasked_intervals.
    """
    results = []
    for m_start, m_end in masked_intervals:
        start = m_start
        for u_start, u_end in unmasked_intervals:
            # no overlap
            if u_end < start or u_start > m_end:
                continue
            # overlap -> cut out overlap
            if u_start > start:
                results.append((start, u_start - 1))
            start = max(start, u_end + 1)
            if start > m_end:
                break
        if start <= m_end:
            results.append((start, m_end))
    return results


def parse_mutations(lines):
    """Return list of (base, pos) or (base, pos, length)."""
    entries = []
    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue
        base, pos = parts[0], int(parts[1])
        if len(parts) == 3:
            length = int(parts[2])
            entries.append((base, pos, length))
        else:
            entries.append((base, pos, False))
    return entries


def generalized_hamming(seq1, seq2, seq_1_name, seq_2_name):
    """
    Compute a generalized Hamming distance between two sequences, considering IUPAC codes.
    Returns the proportion of matches (matches / length).
    """

    IUPAC = {
        "A": {"A"},
        "C": {"C"},
        "G": {"G"},
        "T": {"T"},
        "R": {"A", "G"},
        "Y": {"C", "T"},
        "S": {"G", "C"},
        "W": {"A", "T"},
        "K": {"G", "T"},
        "M": {"A", "C"},
        "B": {"C", "G", "T"},
        "D": {"A", "G", "T"},
        "H": {"A", "C", "T"},
        "V": {"A", "C", "G"},
        "N": {"A", "C", "G", "T", "-"},
        "-": {"-"},
    }

    matches = 0
    mismatches = 0
    i = 0
    for a, b in zip(seq1, seq2):
        i += 1
        if IUPAC[a.upper()] & IUPAC[b.upper()]:  # overlap
            matches += 1
        else:
            mismatches += 1
    if i != 0:
        return (
            matches / i,
            mismatches,
        )  # Return proportion of matches (Goal is to maximize this value) and number of mismatches
    else:
        print(seq_1_name, seq1, len(seq1))
        print(seq_2_name, seq2, len(seq2))
        raise ValueError("Sequences must be of the same non-zero length")
    
def heterozygosity_score(sample_name, cont_seq):
    # Retrieve heterozygous sites for the sample and check if the minor allele is present in the contaminant sequence
    
    # Check that cont_seq is correct
    assert isinstance(cont_seq, list), "cont_seq must be a list"
    # correct length
    assert len(cont_seq) == 29903, "cont_seq must be of length 29903"
    
    
    het_sites = het_sites_dict.get(sample_name, [])
    count = 0
    
    for pos, minor_allele, _ in het_sites:
        # Pos is 1-based in het_sites, convert to 0-based (pos in [1, 29903], cont_seq index in [0, 29902])
        if cont_seq[pos - 1].lower() == minor_allele.lower():
            count += 1

    het_score = count / len(het_sites) if het_sites else 0
    # One in twenty chances to print the het score for debugging
    if random.random() < 0.05:
        # Print sample name and heterozygous positions
        print(f"Sample: {sample_name}, Het score: {het_score:.4f}, Het sites: {len(het_sites)}", flush=True)
    return het_score

def look_for_closest_variant_optimized(
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
            "matching_score",
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
            het_score = 0
            matching_score = -1
            mismatches = -1

            if sample_mut:
                # Expand mutations masked by gen_maple_file
                sample_mut = set(sample_mut.split(";"))
                expanded_mut_set = set()
                for mut in sample_mut:
                    expanded_mut_set.update(expand_ambiguous_mutation(mut, remove_starting_nt=True))
                # Add heterozygous sites to the expanded mutation set
                het_sites = het_sites_dict.get(sample_name, [])
                het_sites_set = set()
                for pos, minor_allele, _ in het_sites:
                    het_mut = f"{pos}{minor_allele}".upper()
                    het_sites_set.update(het_mut)
                expanded_mut_set.update(het_sites_set)
                

                # Find candidate variants
                candidate_counts = Counter()
                for mut in expanded_mut_set:
                    candidate_counts.update(mutation_to_variants.get(mut, []))

                if candidate_counts:
                    max_closeness = max(candidate_counts.values())
                    best_variants = [
                        v
                        for v, count in candidate_counts.items()
                        if count == max_closeness
                    ]

                    # Compute masked regions
                    maple_entry_unmasked = maple_vcf_dict.get(
                        (sample_name, "unmasked"), []
                    )  # UNMASKED ENTRY
                    maple_entry_masked = maple_vcf_dict.get(
                        (sample_name, masked_or_random), []
                    )  # MASKED OR RANDOM ENTRY

                    masked_regions_masked = parse_intervals(
                        maple_entry_masked
                    )  # MASKED REGIONS IN THE MASKED OR RANDOM ENTRY
                    masked_regions_unmasked = parse_intervals(
                        maple_entry_unmasked
                    )  # MASKED REGIONS IN THE UNMASKED (VIRIDIAN) ENTRY
                    masked_regions_alg = subtract_intervals(
                        masked_regions_masked, masked_regions_unmasked
                    )  # REGIONS MASKED BY OUR PROCESS ONLY

                    # Build inverse sequence that contains only the masked regions 
                    viridian_seq = ref_seq.copy()
                    for mut in parse_mutations(maple_entry_unmasked):
                        if len(mut) == 3:
                            base, pos, length = mut
                            viridian_seq[pos - 1 : pos - 1 + length] = [base] * length
                        else:
                            base, pos = mut
                            viridian_seq[pos - 1] = base
                    # THE SEQUENCE IS NOW IDENTICAL TO THE UNMASKED SEQUENCE
                    # Build auxiliary sequence for masked regions
                    aux_viridian_seq = "".join(
                        "".join(viridian_seq[start - 1 : end])
                        for start, end in masked_regions_alg
                    )  # FIRST JOIN TO GET THE CONTIGS THEN JOIN THE CONTIGS

                    # Compare candidates
                    best_contaminant_score = -1
                    # THE SCORE CORRESPONDS TO THE PROPORTION OF MATCHES (matches / length)
                    best_contaminant = ""

                    # Downsample best_variants to at most 100 if too many
                    if len(best_variants) > 100:
                        variants_to_consider = random.sample(best_variants, 100)
                    else:
                        variants_to_consider = best_variants
                    # For each of those (up to 100 candidates), we retrieve their sequence at the masked regions and compute the hamming distance to check how it matches the unmasked sequence
                    for cont_name in variants_to_consider:
                        maple_entry_cont = maple_vcf_dict.get(
                            (cont_name, "clean"), []
                        )  # CONTAMINANT ENTRY IN THE CLEAN TREE
                        cont_seq = ref_seq.copy()
                        for mut in parse_mutations(maple_entry_cont):
                            if len(mut) == 3:
                                base, pos, length = mut
                                cont_seq[pos - 1 : pos - 1 + length] = [base] * length
                            else:
                                base, pos = mut
                                cont_seq[pos - 1] = base
                        

                        aux_cont_seq = "".join(
                            "".join(cont_seq[start - 1 : end])
                            for start, end in masked_regions_alg
                        )

                        hamming_score, n_mismatches = generalized_hamming(
                            aux_viridian_seq, aux_cont_seq, sample_name, cont_name
                        )

                        het_score = heterozygosity_score(sample_name, cont_seq)

                        curr_score = (hamming_score + het_score)/2  # AVERAGE OF HAMMING SCORE AND HET SCORE

                        if curr_score > best_contaminant_score:
                            best_contaminant_score = curr_score
                            best_het_score = het_score
                            min_mismatches = n_mismatches
                            best_contaminant = cont_name
                        # HANDLE TIES
                        elif curr_score == best_contaminant_score:
                            best_contaminant += f";{cont_name}"

                    # Save computed values
                    closest_variant = best_contaminant
                    n_candidates = len(best_variants)
                    closeness = max_closeness
                    closeness_ratio = max_closeness / len(sample_mut)
                    matching_score = best_contaminant_score
                    mismatches = min_mismatches

            # Write updated row
            writer.writerow(
                cols
                + [
                    closest_variant,
                    n_candidates,
                    closeness,
                    closeness_ratio,
                    matching_score,
                ]
            )
    print(f"Updated TSV written to: {output_tsv}", flush=True)

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
            
            masked_positions = {
                pos 
                for start, end in masked_regions_alg 
                for pos in range(start, end + 1)
                }
            
            # Only keep the VCF entries of unmasked_vcf_entry that fall within masked_regions_alg
            alg_masked_mut_set = set()
            alg_masked_mut_set.update(
                f"{pos}{base}".upper()
                for base, pos, is_del_or_n in parse_mutations(unmasked_vcf_entry)
                if not is_del_or_n
                and pos in masked_positions
            )  # if entry is a deletion or Ns, we don't consider it. example: (n, 8765, 14), should we consider deletions?

            # Add heterozygous sites in the entries
            sample_het_sites = het_sites_dict.get(sample_name, [])
            het_mut_set = set()
            het_mut_set.update(
                f"{pos}{minor_allele}".upper() 
                for pos, minor_allele, _ in sample_het_sites
                )

            # alg_masked_vcf_entry now contains only the mutations in the unmasked entry that fall within the regions masked by our process, plus the heterozygous sites
            # This is everything we want to consider for finding the closest variants

            # Expand mutations masked by gen_maple_file
            expanded_pos_set = set()
            expanded_mut_set = set()
            for mut in alg_masked_mut_set:
                expanded_mut = expand_ambiguous_mutation(mut, remove_starting_nt=False)
                if len(expanded_mut) > 1:
                    # We're in the case of an ambiguous mutation, keep this position in memory
                    pos = int(mut[:-1])
                    expanded_pos_set.add(pos)
                expanded_mut_set.update(expanded_mut)
            # Find candidate variants by taking the clean samples that have the most matching entries
            # This method doesn't account for mutations in the contaminants that shouldn't be there

            # STEP 1: Find candidates based on shared mutations
            candidate_counts = Counter()
            candidate_counts.update(chain.from_iterable(mutation_to_variants.get(mut, []) for mut in expanded_mut_set))
            if candidate_counts:
                max_count = candidate_counts.most_common(1)[0][1]
                # max_count should be at least 3
                if max_count >= 3:
                    
                    best_candidates = [
                        v
                        for v, c in candidate_counts.items()
                        if c == max_count
                    ]

                    # STEP 2: Iterate on the candidates, check the matching between their vcf entries and the alg_masked_vcf_entry on the masked regions and on the heterozygous sites
                    for candidate in best_candidates:
                        # Retrieve the vcf entry of the candidate
                        candidate_vcf_entry = maple_vcf_dict.get(
                            (candidate, "clean"), []
                        )
                        candidate_mut_set = {
                            f"{pos}{base}".upper()  # Format it to a string like 2343T
                            for base, pos, is_del_or_n in parse_mutations(candidate_vcf_entry)
                            if not is_del_or_n
                        }  # If entry is a deletion or Ns, we don't consider it. example: (n, 8765, 14)
                        candidate_mut_on_masked_regions_set = {
                            m
                            for m in candidate_mut_set
                            if int(m[:-1]) in masked_positions
                        }  # Only add the mutation if it falls within the masked regions

                        # Check intersection, penalize missing mutations and extra mutations
                        shared = expanded_mut_set & candidate_mut_on_masked_regions_set
                        extra = candidate_mut_on_masked_regions_set - expanded_mut_set
                        shared_het = het_mut_set & candidate_mut_set

                        score_1 = (len(shared) - len(extra)) / (len(shared) + len(extra) or 1) # Ranges in [-1, 1]
                        score_1 = (score_1 + 1) / 2  # Ranges in [0, 1]
                        score_2 = len(shared_het) / (len(het_mut_set) or 1) # Ranges in [0, 1]

                        ratio_score = (2*score_1 + score_2) / 3

                        if ratio_score > closeness_ratio:
                            best_contaminant = candidate
                            closeness_ratio = ratio_score
                            closeness = len(shared)
                            # generate a hash of the shared and shared_het sets to avoid duplicates
                            best_hash = [hash((frozenset(shared), frozenset(shared_het)))]

                        elif ratio_score == closeness_ratio:
                            # check that current candidate doesn't have exactly the same shared and shared het mutations as previous best candidates
                            candidate_hash = hash((frozenset(shared), frozenset(shared_het)))
                            # Values before hash Min n_candidates: 0 Max n_candidates: 115326 Mean n_candidates: 229.8029776674938
                            if candidate_hash not in best_hash:
                                best_contaminant += f";{candidate}"
                                best_hash.append(candidate_hash)
                    # Save computed values
                    closest_variant = best_contaminant
                    n_candidates = len(best_contaminant.split(";"))
                
                else:
                    closest_variant = "Not enough shared mutations with any candidate"
                    n_candidates = 0
                    closeness = 0
                    closeness_ratio = 0

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

def parallelize_contaminant_search(
    input_tsv,
    output_tsv,
    mutation_col_name,
    masked_or_random,
    ref_seq,
    num_cores,
    multiprocess=True,
):
    """Parallelize the contaminant search across multiple cores.
    Split the input TSV into chunks and process each chunk in parallel.
    """
    
    if multiprocess:
    
        max_cores = cpu_count()
        if num_cores > max_cores:
            print(f"Requested {num_cores} cores, but only {max_cores} available. Using {max_cores}.", flush=True)
            num_cores = max_cores
        print(f"Using {num_cores} cores for parallel processing.", flush=True)
        # Split the input TSV into chunks
        input_chunk_files = []
        output_chunk_files = []
        with open(input_tsv, "r") as f_in:
            header = f_in.readline() # Read header
            lines = f_in.readlines()
            chunk_size = len(lines) // num_cores + 1
            for i in range(num_cores):
                chunk_file_path = Path(str(input_tsv).replace(".tsv", f"_chunk_{i}.tsv"))
                input_chunk_files.append(chunk_file_path)
                output_chunk_file_path = Path(str(output_tsv).replace(".tsv", f"_chunk_{i}.tsv"))
                output_chunk_files.append(output_chunk_file_path)
                with open(chunk_file_path, "w") as f_chunk:
                    f_chunk.write(header)
                    f_chunk.writelines(lines[i*chunk_size:(i+1)*chunk_size])
            print(f"Created {len(input_chunk_files)} chunk files.", flush=True)

        # Process each chunk in parallel using look_for_closest_variant_optimized
        # Prepare arguments for each process
        args = []
        for i in range(num_cores):
            args.append((
                str(input_chunk_files[i]),
                str(output_chunk_files[i]),
                mutation_col_name,
                masked_or_random,
                ref_seq,
            ))

        with Pool(processes=num_cores) as pool:
            pool.starmap(look_for_closest_variant_optimized_unified, args)
        
        print("All chunks processed, combining results...", flush=True)
        # Combine chunk outputs into final output TSV
        with open(output_tsv, "w") as f_out:
            f_out.write(header)
            for output_chunk in output_chunk_files:
                with open(output_chunk, "r") as f_chunk:
                    next(f_chunk)  # Skip header
                    f_out.writelines(f_chunk.readlines())
                    
                # Remove chunk files
                output_chunk.unlink()

        for input_chunk in input_chunk_files:
            input_chunk.unlink()
            
    else:
        # Single process
        look_for_closest_variant_optimized_unified(
            input_tsv,
            output_tsv,
            mutation_col_name,
            masked_or_random,
            ref_seq,
        )

    print(f"Final output written to: {output_tsv}", flush=True)


def load_globals():
    """Load global dictionaries from pickle files or build them if they don't exist.
    """

    global mutation_to_variants, maple_vcf_dict, het_sites_dict

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
        
    with open(het_sites_dict_path, "rb") as f:
        het_sites_dict = pickle.load(f)
    print(f"Loaded heterozygous sites dictionary with {len(het_sites_dict)} samples.", flush=True)


if __name__ == "__main__":

    load_globals()

    # GET REF SEQUENCE
    with open(path_ref_seq, "r") as ref_file:
        ref_file.readline()  # Skip header
        ref_seq = ref_file.readline().strip()
        ref_seq = list(ref_seq)

    if masked_or_random == "masked":
        path_store = path_psmp_masked

    elif masked_or_random == "random":
        path_store = path_psmp_random

    print(f"Processing {masked_or_random} samples...")
    path_store_updated = Path(
        str(path_store).replace("results", "results_with_contaminants")
    )
    
    # Add contaminant candidates
    parallelize_contaminant_search(
        input_tsv=path_store,
        output_tsv=path_store_updated,
        mutation_col_name="mutations_masked",
        masked_or_random=masked_or_random,
        ref_seq=ref_seq,
        num_cores=num_cores,
        multiprocess=False,
    )

    params_path = "/nfs/research/goldman/anoufa/src/dpca/_params.py"

    update_params_file(
        params_path,
        {
            f"processed_placement_{masked_or_random}_with_contaminants": str(
                path_store_updated
            )
        },
    )
    print(f"Processed placements results saved to {path_store_updated}")
