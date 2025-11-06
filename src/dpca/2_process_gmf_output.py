from pathlib import Path
import pandas as pd
from tqdm import tqdm
import random
import shutil
import gzip
import _params
from _aux_functions import generate_sh_param_file, compress_file, build_maple_entry, build_maple_file, update_params_file, expand_ambiguous_mutation
from multiprocessing import Pool, cpu_count
import os
from functools import partial
import pickle

data_dir = _params.data_dir
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
param_term = _params.param_term
num_cores = _params.num_cores
masking_method = _params.masking_method




def build_mut_dict(maple_file):
    """Build a dictionary with keys sample names and values mutations.
    The mutations are written as strings in the format "100T", meaning position 100 has nucleotide T.
    The nucleotide present before the mutation is removed.
    The dictionary has the following structure:
    {
        "sample_name_1": {"100T", "200G", ...},
        "sample_name_2": {"150A", "300C", ...},

    Args:
        maple_file (str): Path to the maple file.

    Returns:
        dict: A dictionary with sample names as keys and sets of mutations as values.
    """
    variant_mut_dict = {}
    seen_mutation_tuples = set() # To store unique sorted mutation tuples
    
    with open(maple_file, "r") as f:
        # Read the first line (header)
        f.readline()
        # Read the second line (reference sequence)
        ref_seq_line = f.readline().strip()
        ref_seq = ref_seq_line 
        
        current_variant_mutations = []
        seq_name = "" 

        # Using tqdm with f for progress bar for lines
        for line in tqdm(f, desc="Building mutation dictionary", mininterval=60):
            line = line.strip()

            if line.startswith(">"):
                if current_variant_mutations: # Check if the previous variant had mutations
                    # Process the previous variant's mutations
                    # Sort the mutations and convert to a set for hashing
                    current_variant_mutations = set(current_variant_mutations)
                    expanded_mutations = set()
                    for mut in current_variant_mutations:
                        expanded_mutations.update(expand_ambiguous_mutation(mut))

                    variant_mut_dict[seq_name] = expanded_mutations

                # Start a new sequence
                seq_name = line # The full line starting with '>' is the sequence name
                # Remove the '>' character
                seq_name = seq_name[1:]
                current_variant_mutations = [] # Reset for the new variant
            else:
                # Example line format: "g\t 3456" or "-\t 7890\t 3"            
                # Check if the line looks like a mutation line
                parts = line.split('\t')
                if parts[0].upper() in ["A", "C", "G", "T"]:
                    new_nt, pos = parts[0], parts[1] 
                    mut = str(pos) + str(new_nt)
                    mut = mut.upper()  # Ensure the mutation is in uppercase
                    current_variant_mutations.append(mut)

        # After the loop, save the very last sequence
        if current_variant_mutations:
            current_variant_mutations = set(current_variant_mutations)
            expanded_mutations = set()
            for mut in current_variant_mutations:
                expanded_mutations.update(expand_ambiguous_mutation(mut))

            variant_mut_dict[seq_name] = expanded_mutations
                
    return variant_mut_dict

def append_clean_tree(align_file, content):
    # print(f"[Process {os.getpid()}] Appending clean samples to {align_file}")
    with open(align_file, "ab") as f:
        f.write(content)

if __name__ == "__main__":
    # Generate the batchs
                                                
    # Access the total median cat list and add the results to the file
    final_mut_path = Path(f"{data_dir}/2/n_masked_and_masked_mut/masked_mut_{param_term}.tsv")

    final_mut_tsv = pd.DataFrame(columns=["sample_name", "masked_mutations_masked", "masked_mutations_random"])

    batches_folder_path = Path(f"{data_dir}/1/")

    for mut_tsv_file in batches_folder_path.glob(f"maple_mutations_batch*_{param_term}.tsv"):

        # Load the data
        mut_tsv = pd.read_csv(mut_tsv_file, sep="\t")
    
        # Add the data to the final dataframe
        final_mut_tsv = pd.concat([final_mut_tsv, mut_tsv], ignore_index=True)
        mut_tsv_file.unlink()
        
    print(f"Processed {len(final_mut_tsv)} mutation entries.")
    
    final_mut_tsv.to_csv(final_mut_path, sep="\t", index=False)
    

    final_clean_tree_path = f"{data_dir}/2/alignment_files/clean_tree_alignment_file_{param_term}.maple"
    maple_alignment_unmasked_file_path = f"{data_dir}/2/alignment_files/unmasked_alignment_{param_term}.maple"
    maple_alignment_random_file_path = f"{data_dir}/2/alignment_files/random_alignment_{param_term}.maple"
    maple_alignment_masked_file_path = f"{data_dir}/2/alignment_files/masked_alignment_{param_term}.maple"
    build_maple_file(path_ref_seq, maple_alignment_unmasked_file_path)
    build_maple_file(path_ref_seq, maple_alignment_random_file_path)
    build_maple_file(path_ref_seq, maple_alignment_masked_file_path)
    build_maple_file(path_ref_seq, final_clean_tree_path)
    
    clean_samples_maple_list = []
    current_sample_entry = []

    for clean_tree_file in batches_folder_path.glob(f"clean_tree_batch*_{param_term}.maple"):
        with open(clean_tree_file, "r") as f:
            # Add each entry as a list to the clean_samples_maple_list
            for line in f:
                if line.startswith(">"):
                    if current_sample_entry: # If not the first entry, save the previous one
                        clean_samples_maple_list.append("\n".join(current_sample_entry) + "\n")
                    current_sample_entry = [line.strip()] # Start a new entry
                else:
                    current_sample_entry.append(line.strip()) # Add to the current entry
            if current_sample_entry: # Save the last entry
                clean_samples_maple_list.append("\n".join(current_sample_entry) + "\n")
                current_sample_entry = []
                
        clean_tree_file.unlink()

    # Remove identical entries (if any)
    print(f"Number of clean samples before removing duplicates: {len(clean_samples_maple_list)}")
    unique_samples_list = list(set(clean_samples_maple_list))
    
    with open(final_clean_tree_path, "a") as new_f:
        new_f.writelines(unique_samples_list)
    print(f"Number of clean samples after removing duplicates: {len(unique_samples_list)}")
    print(f"Clean tree file saved")

    for cons_align_file in batches_folder_path.glob(f"unmasked_alignment_batch*_{param_term}.maple"):
        # Copy the file content to the MAPLE input folder
        with open(cons_align_file, "r") as f:
            with open(maple_alignment_unmasked_file_path, "a") as new_f:
                new_f.write(f.read())
                
        cons_align_file.unlink()


    for random_align_file in batches_folder_path.glob(f"random_alignment_batch*_{param_term}.maple"):
        # Copy the file content to the MAPLE input folder
        with open(random_align_file, "r") as f:
            with open(maple_alignment_random_file_path, "a") as new_f:
                new_f.write(f.read())
                
        random_align_file.unlink()
        

    for masked_align_file in batches_folder_path.glob(f"masked_alignment_batch*_{param_term}.maple"):
        # Copy the file content to the MAPLE input folder
        with open(masked_align_file, "r") as f:
            with open(maple_alignment_masked_file_path, "a") as new_f:
                new_f.write(f.read())
                
        masked_align_file.unlink()

    print(f"Alignment files saved")

    final_path_n_masked = Path(f"{data_dir}/2/n_masked_and_masked_mut/n_masked_{param_term}.tsv")
    
    # Open the final output file in write mode.
    # This will create the file if it doesn't exist or overwrite it if it does.
    with open(final_path_n_masked, 'w') as outfile:
        # Iterate over all .tsv files in the specified folder
        for n_masked_tsv_file in tqdm(batches_folder_path.glob(f"n_masked_batch_*_{param_term}.tsv"), mininterval=30, desc="Processing n_masked files"):
            with open(n_masked_tsv_file, 'r') as infile:
                lines = infile.readlines()

                outfile.writelines(lines)

            n_masked_tsv_file.unlink()

    print(f"n_masked files saved")
    
    # Count number of samples in clean tree file
    n_samples = 0
    with open(final_clean_tree_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                n_samples += 1
    print(f"Number of samples in clean tree file: {n_samples}")

    # CHECK IF THERE IS AN ALIGNMENT FILE TO APPEND IN THE CLEAN TREE FOLDER, OTHERWISE KEEP THE ONE GENERATED HERE
    potential_existing_alignment_file = Path(f"{data_dir}/2/clean_tree/clean_tree_alignment_file_{param_term}.maple")
    if potential_existing_alignment_file.exists():
        print(f"An existing alignment file was found at {potential_existing_alignment_file}. It will be used instead of the newly generated one.")
        final_clean_tree_path = potential_existing_alignment_file
        
        
    # Generate the alignment files to prepare the sample placement
    with open(final_clean_tree_path, "r") as f:
        next(f)
        next(f) # Skip the first two lines (reference sequence)
        clean_tree_content = f.read().encode()

    align_files = list(batches_folder_path.glob(f"maple_alignment_{masking_method}_batch*_{param_term}.maple"))
    print(f"Appending clean samples to {len(align_files)} alignment files using {num_cores} CPU cores...")

    with Pool(processes=num_cores) as pool:
        pool.map(partial(append_clean_tree, content=clean_tree_content), align_files)
            
    print("Added clean samples to alignment files")
    
    # Build mutation dictionary
    var_mut_dict = build_mut_dict(final_clean_tree_path) # ref_seq is read from the file
    print(f"Built mutation dictionary with {len(var_mut_dict)} unique variant entries.")
    
    # Build mut_to_variant_dict, where keys are mutations and values are lists of variants having that mutation
    mut_to_variant_dict = {}
    for variant, mutations in var_mut_dict.items():
        for mut in mutations:
            if mut not in mut_to_variant_dict:
                mut_to_variant_dict[mut] = []
            mut_to_variant_dict[mut].append(variant)
            
    # sanity check, print the first entries of both dictionaries
    print("Sample entries from variant to mutations dictionary:")
    for i, (variant, mutations) in enumerate(var_mut_dict.items()):
        # mutation is a set, show some elements
        print(f"{variant}: {list(mutations)[:5]}")  # Print only first 5 mutations for brevity
        if i >= 2:
            break
    print("\nSample entries from mutation to variants dictionary:")
    for i, (mut, variants) in enumerate(mut_to_variant_dict.items()):
        print(f"{mut}: {list(variants)[:5]}")  # Print only first 5 variants for brevity
        if i >= 2:
            break        
    
    # Save the variant:mutation dictionary to a pickle file
    dict_path = f"{data_dir}/2/dict/clean_samples_var_to_mut_dict_{param_term}.pickle"
    with open(dict_path, "wb") as f:
        pickle.dump(var_mut_dict, f)
    print(f"Variant to mutation dictionary saved to {dict_path}.")
    # Save the mutation:variant dictionary to a pickle file
    inv_dict_path = f"{data_dir}/2/dict/clean_samples_mut_to_var_dict_{param_term}.pickle"
    with open(inv_dict_path, "wb") as f:
        pickle.dump(mut_to_variant_dict, f)
    print(f"Mutation to variant dictionary saved to {inv_dict_path}.")
    