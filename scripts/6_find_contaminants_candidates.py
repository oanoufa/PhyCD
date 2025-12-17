# This script should go over the tsv results of MAPLE sample placement and
# generate the distributions of the branch lengths for the three types of sequences
# (viridian, masked, random).

# The goal is to see whether the masked samples are closer to the tree than the viridian ones or not.

from pathlib import Path
from tqdm import tqdm
from collections import Counter
import argparse
import csv
import random
from _aux_functions import update_params_file, expand_ambiguous_mutation, generate_sample_list, parse_maple_file, smart_open, save_pickle_dict, load_pickle_dict
from multiprocessing import Pool
import os
from itertools import chain

parser = argparse.ArgumentParser(
    description="Process the concatenated output of the MAPLE sample placement."
)

parser.add_argument("--masked_or_random", type=str,
                    help="Either masked or random, to process both in parallel")
parser.add_argument('--data_dir', type=str,
                    help='Directory containing all the files generated during the pipeline.')
parser.add_argument('--samples_dir', type=str,
                    help='Directory containing the sample names and paths to their Viridian output folders.')
parser.add_argument('--param_term', type=str,
                    help='str built from the initial parameters used to recognize the files.')
parser.add_argument('--path_ref_seq', type=str,
                    help='Path to the reference sequence used for the pipeline.')
parser.add_argument('--num_cores', type=int, default=8,
                    help='Number of cores to use to multiprocess the contaminant search.')
parser.add_argument('--compress', type=int,
                    help='1 or 0, enables compression.')
args = parser.parse_args()

# Get the arguments
masked_or_random = args.masked_or_random
data_dir = args.data_dir
samples_dir = args.samples_dir
param_term = args.param_term
path_ref_seq = args.path_ref_seq
num_cores = args.num_cores
compress = args.compress

# Global variables
mutation_to_variants = None
maple_vcf_dict = None
het_sites_dict = None

def build_maple_vcf_dict(clean_file_path, param_term, masked_or_random):
    """Build a dictionary mapping sample names and type (clean, random, masked, unmasked) to their MAPLE VCF entries."""
    maple_vcf_dict = {}
    # FIRST RETRIEVE ALL SAMPLES FROM CLEAN TREE
    print("Building MAPLE VCF dictionary from clean tree alignment...", flush=True)
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

    alignment_folder = Path(f"{data_dir}/2/alignment_files/")
    print(f"Adding unmasked and {masked_or_random} entries...", flush=True)
    for alignment_file in tqdm(
        alignment_folder.glob(f"*_alignment_{param_term}*"),
        desc="Building MAPLE VCF dictionary",
        mininterval=60,
    ):
        # OPEN THE 3 FILES unmasked_alignment, masked_alignment, random_alignment
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

        with smart_open(alignment_file, "rt") as file:
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

def look_for_best_contaminant(
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
    - best_contaminant: the name(s) of the best contaminant(s) in the clean tree
    - n_candidates: number of candidate variants with the highest closeness score
    - score 1 and 2: closeness scores
    - final_score: average between scores 1 and 2
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
            "best_contaminant",
            "n_candidates",
            "n_unique_candidates",
            "score_1",
            "score_2",
            "final_score",
        ]

        writer.writerow(new_header)

        for cols in tqdm(reader, desc=f"[PID {os.getpid()}] Processing samples", mininterval=120):
            sample_name = cols[sample_idx]
            sample_mut = cols[mutation_idx]

            # Default output values
            best_contaminant = "No mutations masked by gen_maple_file"
            n_candidates = 0
            n_unique_candidates = 0
            best_score_1 = 0
            best_score_2 = 0
            best_final_score = 0
            best_hash = []
            
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
            alg_masked_mut_set = set() # unexpanded mutations of the unmasked sequence
            alg_masked_mut_set.update(
                f"{pos}{base}".upper()
                for base, pos, is_del_or_n in parse_mutations(unmasked_vcf_entry)
                if not is_del_or_n
                and pos in masked_positions
            )  # if entry is a deletion or Ns, we don't consider it. example: (n, 8765, 14), should we consider deletions?
            #alg_masked_mut_set contains mutations in the unmasked entry that fall within the masked regions
            
            # Add heterozygous sites in the entries
            sample_het_sites = het_sites_dict.get(sample_name, [])
            het_mut_set = set()
            het_mut_set.update(
                f"{pos}{minor_allele}".upper() 
                for pos, minor_allele, _ in sample_het_sites
                )
            # het_mut_set contains minor alleles at heterozygous sites

            # Expand mutations masked by gen_maple_file
            expanded_mut_set = set()
            for mut in alg_masked_mut_set:
                expanded_mut = expand_ambiguous_mutation(mut, remove_starting_nt=False)
                expanded_mut_set.update(expanded_mut)
            # Find candidate variants by taking the clean samples that have the most matching entries
            # This method doesn't account for mutations in the contaminants that shouldn't be there

            # STEP 1: Find candidates based on shared mutations
            candidate_counts = Counter()
            candidate_counts.update(chain.from_iterable(mutation_to_variants.get(mut, []) for mut in expanded_mut_set))

            if not candidate_counts: # If no mutations were masked by gen_maple_file, try to match with heterozygous sites
                candidate_counts = Counter()
                candidate_counts.update(chain.from_iterable(mutation_to_variants.get(mut, []) for mut in het_mut_set))
            
            if candidate_counts:
                max_count = candidate_counts.most_common(1)[0][1]
                # max_count should be at least 3
                if max_count >= 0:
                    
                    best_candidates = [
                        v
                        for v, c in candidate_counts.items()
                        if c == max_count
                    ]

                    # STEP 2: Iterate on the candidates and compute scores depending on shared mutations, penalizing extra mutations
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

                        score_1 = len(shared) / (len(alg_masked_mut_set) + len(extra) or 1) # Ranges in [0, 1]
                        score_2 = len(shared_het) / (len(het_mut_set) or 1) # Ranges in [0, 1]

                        if len(het_mut_set) == 0:
                            print(sample_name)
                        
                        final_score = (score_1 + score_2) / 2

                        if final_score > best_final_score:
                            best_contaminant = candidate
                            best_score_1 = score_1
                            best_score_2 = score_2
                            best_final_score = final_score
                            n_candidates = 1
                            n_unique_candidates = 1
                            # generate a hash of the shared and shared_het sets to avoid duplicates
                            best_hash = [hash((frozenset(shared), frozenset(shared_het)))]

                        elif final_score == best_final_score:
                            # check that current candidate doesn't have exactly the same shared and shared het mutations as previous best candidates
                            candidate_hash = hash((frozenset(shared), frozenset(shared_het)))
                            # Values before hash Min n_candidates: 0 Max n_candidates: 115326 Mean n_candidates: 229.8029776674938
                            if candidate_hash not in best_hash:
                                best_contaminant += f";{candidate}"
                                best_hash.append(candidate_hash)
                                n_unique_candidates += 1
                            n_candidates += 1
                
                else:
                    best_contaminant = "Not enough shared mutations with any candidate"
                    n_candidates = 0
                    n_unique_candidates = 0
                    best_score_1 = 0
                    best_score_2 = 0
                    best_final_score = 0

            # Write updated row
            writer.writerow(
                cols
                + [
                    best_contaminant,
                    n_candidates,
                    n_unique_candidates,
                    best_score_1,
                    best_score_2,
                    best_final_score,
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
    
        max_cores = os.cpu_count()
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
            pool.starmap(look_for_best_contaminant, args)
        
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
        look_for_best_contaminant(
            input_tsv,
            output_tsv,
            mutation_col_name,
            masked_or_random,
            ref_seq,
        )

    print(f"Final output written to: {output_tsv}", flush=True)

def prepare_dataset_eyre_model(path_df_with_contaminants, ref_seq, masked_or_random):
    """Generate the files needed for the Eyre et al. model based on the identified contaminants.

    Args:
        path_df_with_contaminants (str): Path to the previously generated dataset
        ref_seq (list): Reference sequence as a list of characters
        masked_or_random (str): Indicates whether the data is masked or random
    """
    output_folder = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/")
    output_folder.mkdir(parents=True, exist_ok=True)
    
    base_counts_folder = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/base_counts/")
    base_counts_folder.mkdir(parents=True, exist_ok=True)
    
    pairs_dict = {}
    pairs_dict_path = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/sample_contaminant_pairs_{masked_or_random}.pickle")

    variable_sites_per_sample_dict = {}
    variable_sites_per_sample_dict_path = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/variable_sites_per_sample_{masked_or_random}.pickle")
    
    # pairs_dict stores sample_name -> set of contaminant names
    maple_entries_dict = {}
    # maple_entries_dict stores sequence_name -> list of lines from the MAPLE file
    with open(path_df_with_contaminants, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in tqdm(reader, desc="Iterating over potentially contaminated samples", mininterval=90):
            sample_name = row["sample_name"]
            contaminant_name = row["best_contaminant"]
            
            # if best_contaminant contains "mutations" it means that no contaminants were found
            if "mutations" in contaminant_name:
                continue
            pairs_dict[sample_name] = set(contaminant_name.split(";"))
            # Add the sample itself to the its set of contaminants as a null hypothesis
            pairs_dict[sample_name].add(sample_name)
            # Add the maple entries of the sample
            maple_entries_dict[sample_name] = maple_vcf_dict.get((sample_name, masked_or_random), [])
            for cont in contaminant_name.split(";"):
                if cont not in maple_entries_dict:
                    maple_entries_dict[cont] = maple_vcf_dict.get((cont, "clean"), [])

            # Store the variable sites per sample
            variable_sites = set()
            for sample_or_cont in pairs_dict[sample_name]:
                entry_lines = maple_entries_dict.get(sample_or_cont, [])
                for line in entry_lines:
                    parts = line.split()
                    pos = int(parts[1])
                    variable_sites.add(pos)
            variable_sites_per_sample_dict[sample_name] = variable_sites
    
    # Store pairs_dict for the eyre model
    save_pickle_dict(pairs_dict, pairs_dict_path, compress)
    
    # Reconstruct full sequences
    reconstructed_seqs = {}
    for seq_name, entry_lines in tqdm(maple_entries_dict.items(), desc="Reconstructing sequences at variable sites", mininterval=90):
        seq = list(ref_seq)
        for line in entry_lines:
            parts = line.split()
            pos = int(parts[1])
            base = parts[0].upper()
            if base in ["A", "C", "G", "T"]:
                seq[pos - 1] = base  # MAPLE positions are 1-based
        
        # full sequences
        reconstructed_seqs[seq_name] = ("".join(seq)).upper()
    
    # Write the reconstructed sequences to the haplotype_sequence.txt file, with two columns: sequence_name, sequence (unnamed)
    haplotype_sequence_path = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/haplotype_sequence.txt")
    with haplotype_sequence_path.open("w") as f:
        for seq_name, seq in reconstructed_seqs.items():
            f.write(f"{seq_name}\t{seq}\n")
            
    # Generate the haplotype_frequency.txt file, with two columns: sequence_name, frequency
    # For simplicity, we will assign a frequency of 1/number_of_sequences to each sequence
    haplotype_frequency_path = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/haplotype_frequency.txt")
    frequency = 1.0
    with haplotype_frequency_path.open("w") as f:
        for seq_name in reconstructed_seqs.keys():
            f.write(f"{seq_name}\t{frequency}\n")
        
    # Clear previous samples
    for file in base_counts_folder.glob("*_base_counts.txt"):
        file.unlink()

    # Go look for the qc files of the samples and write the counts of nt at variable sites only of those samples
    potentially_contaminated_samples_set = set(pairs_dict.keys())
    name_path_list = generate_sample_list(potentially_contaminated_samples_set, samples_dir=samples_dir)
    
    for path, sample_name in tqdm(name_path_list, desc="Processing samples to write base counts", unit="sample", mininterval=120):
        
        ref_pos_list = []
        unma_seq = []
        counts_list = []
        
        variable_sites = variable_sites_per_sample_dict.get(sample_name, set())
        
        with smart_open(path, "rt") as f:  # read as text directly
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
                    
                    n_0 = sum(1 for count in counts if count == 0)
                    # list like ref_pos, A_count, C_count, G_count, T_count
                    if ref_pos in variable_sites:
                    # check that at least one count > 0
                        if sum(counts) == 0:
                            # remove this site from the variable sites of the sample
                            variable_sites.remove(ref_pos)
                        else:
                            counts_list.append([ref_pos] + counts)
                    
                else:
                    # Deletion ( - row) or rows full of dots (start/end gaps)
                    ref_pos = int(line[0])
                    cons_nt = line[4].lower()
                    if ref_pos in variable_sites:
                        # remove site from variable sites
                        variable_sites.remove(ref_pos)
                    
                # track unmasked sequence
                if ref_pos_list and ref_pos == ref_pos_list[-1]:
                    unma_seq[-1] += cons_nt
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
        
        # Save the new version of the variable sites of this sample
        variable_sites_per_sample_dict[sample_name] = variable_sites
                    
        # Now write the counts at variable sites only
        output_counts_path = base_counts_folder / f"{sample_name}_base_counts.txt"
        with output_counts_path.open("w") as f:
            # header
            f.write("site\tA\tC\tG\tT\n")
            for counts in counts_list:
                f.write("\t".join(map(str, counts)) + "\n")

    # Write the list of file names to filenames_list.txt
    filenames_list_path = Path(f"{data_dir}/6/eyre_model_{masked_or_random}/filenames_list.txt")
    with filenames_list_path.open("w") as f:
        for file_name in base_counts_folder.glob("*_base_counts.txt"):
            f.write(f"base_counts/{file_name.name}\n")
    
    # Store variable_sites_per_sample_dict for the eyre model
    save_pickle_dict(variable_sites_per_sample_dict, variable_sites_per_sample_dict_path, compress)

    # Variable sites are defined per sample and are defined as:
    # the sites where the sample or one of its potential contaminant have a maple entry
    # AND sites where the sample has count data (not 0 everywhere)

def load_globals():
    """Load global dictionaries from tsv files or build them if they don't exist.
    """

    global mutation_to_variants, maple_vcf_dict, het_sites_dict

    # Load the mutation to variants dictionary once
    mutation_to_variants_path = Path(
        f"{data_dir}/2/dict/clean_samples_mut_to_var_dict_{param_term}.pickle"
    )

    clean_tree_alignment_path = Path(
        f"{data_dir}/2/alignment_files/clean_tree_alignment_file_{param_term}.maple"
    )

    maple_vcf_dict_path = Path(
        f"{data_dir}/2/dict/maple_vcf_dict_{param_term}_{masked_or_random}.pickle"
    )

    het_sites_dict_path = Path(f"{data_dir}/2/dict/het_sites_dict_{masked_or_random}_{param_term}.pickle")

    # BUILD MAPLE VCF DICT
    maple_vcf_dict = build_maple_vcf_dict(clean_tree_alignment_path, param_term, masked_or_random)
    print("Storing MAPLE VCF dictionary...", flush=True)
    save_pickle_dict(maple_vcf_dict, maple_vcf_dict_path, compress)
    print(f"Stored MAPLE VCF dictionary with {len(maple_vcf_dict)} entries to {maple_vcf_dict_path}.", flush=True)

    # LOAD MUT TO VAR AND HET SITES DICT
    mutation_to_variants = load_pickle_dict(mutation_to_variants_path, compress)
    print(f"Loaded existing mutation to variants dictionary with {len(mutation_to_variants)} mutations.", flush=True)

    het_sites_dict = load_pickle_dict(het_sites_dict_path, compress)
    print(f"Loaded heterozygous sites dictionary with {len(het_sites_dict)} samples.", flush=True)


if __name__ == "__main__":

    load_globals()

    # GET REF SEQUENCE
    with open(path_ref_seq, "r") as ref_file:
        ref_file.readline()  # Skip header
        ref_seq = ref_file.readline().strip()
        ref_seq = list(ref_seq)
        
    path_psmp = Path(f"{data_dir}/5/processed_placements_results_{masked_or_random}_{param_term}.tsv")

    print(f"Processing {masked_or_random} samples...")
    path_psmp_updated = Path(
        (str(path_psmp).replace("results", "results_with_contaminants")).replace("/5/", "/6/")
    )
    
    # Add contaminant candidates
    parallelize_contaminant_search(
        input_tsv=path_psmp,
        output_tsv=path_psmp_updated,
        mutation_col_name="mutations_masked",
        masked_or_random=masked_or_random,
        ref_seq=ref_seq,
        num_cores=num_cores,
        multiprocess=False,
    )
    
    # Generate the data for the adapted_eyre_model
    prepare_dataset_eyre_model(
        path_df_with_contaminants=path_psmp_updated,
        ref_seq=ref_seq,
        masked_or_random=masked_or_random,
    )
