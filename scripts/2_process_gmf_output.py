from pathlib import Path
import pandas as pd
from tqdm import tqdm
import shutil
from _aux_functions import compress_file, build_maple_entry, build_maple_file, expand_ambiguous_mutation, smart_open, save_pickle_dict
from multiprocessing import Pool, cpu_count
from functools import partial
import argparse
import numpy as np
import plotly.express as px
import zstandard as zstd

parser = argparse.ArgumentParser(description='Process the output of the first batched script.')
parser.add_argument('--data_dir', type=str,
                    help='Directory containing all the files generated during the pipeline.')
parser.add_argument('--param_term', type=str,
                    help='str built from the initial parameters used to recognize the files.')
parser.add_argument('--path_ref_seq', type=str,
                    help='Path to the reference sequence used for the pipeline.')
parser.add_argument('--path_metadata_tsv', type=str,
                    help='Path to the TSV of Viridian samples metadata.')
parser.add_argument('--num_cores', type=int, default=8,
                    help='Number of cores to use to append the clean samples to the alignment files.')
parser.add_argument('--compress', type=int,
                    help='1 or 0, enables compression.')
args = parser.parse_args()

# Get the arguments
data_dir = args.data_dir
param_term = args.param_term
path_ref_seq = args.path_ref_seq
num_cores = args.num_cores
path_metadata_tsv = args.path_metadata_tsv
compress = args.compress


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
                # Example line format: "g\t3456" or "-\t7890\t3"            
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
    """Append `content` to align_file only if it is not already the last block."""
    align_file = Path(align_file)
    if isinstance(content, str):
        content = content.encode()

    # Instead of seeking, read sequentially and check if content is in last N bytes
    N = len(content) + 1000  # a buffer in case content is not exactly at the end
    existing_end = b""
    with smart_open(align_file, "rb") as f:
        while chunk := f.read(65536):
            existing_end = (existing_end + chunk)[-N:]

    if content in existing_end:
        print(f"{align_file} already contains content — skipping.")
        return

    # Append
    if align_file.suffix == ".zst":
        cctx = zstd.ZstdCompressor()
        with open(align_file, "ab") as f:  # raw file
            with cctx.stream_writer(f) as writer:
                writer.write(content)
    else:
        with smart_open(align_file, "ab") as f:
            f.write(content)

def plot_reference_samples_stacked_barplot(n_masked_path_tsv, data_dir, param_term):
    """Plot the stacked bar plot showing the samples binned by dropout masked positions, stacked by the number of heterozygous sites in the sample.
    """
    df_n_masked = pd.read_csv(n_masked_path_tsv, sep='\t')
    df_n_masked['dropout_masked'] = df_n_masked['n_masked_masked'] - df_n_masked['n_masked_cons']

    # Clip values
    max_dropout = 1000
    n_bins = 20
    step = int(max_dropout/n_bins)
    max_het = 8
    dropout_bins = [0, 1, 200, 400, 600, 1000, 1500, 2000, np.inf]
    labels_dropout = ['0'] + [f"{dropout_bins[i]} - {dropout_bins[i+1]-1}" for i in range(1, len(dropout_bins)-2)] + [f"{dropout_bins[-2]}+"]
    het_bins = [i for i in range(max_het+1)] + [np.inf]
    labels_het = [str(i) for i in range(max_het)] + [f"{max_het}+"]

    df_n_masked['dropout_bin'] = pd.cut(df_n_masked['dropout_masked'], bins=dropout_bins, labels=labels_dropout, right=False)
    df_n_masked['het_bin'] = pd.cut(df_n_masked['n_het_sites'], bins=het_bins, labels=labels_het, right=False)
    df_n_masked['dropout_bin'] = df_n_masked['dropout_bin'].astype(str)
    df_n_masked['het_bin'] = df_n_masked['het_bin'].astype(str)

    df_n_masked['dropout_bin'] = pd.Categorical(
        df_n_masked['dropout_bin'],
        categories=labels_dropout,
        ordered=True
    )

    ct = pd.crosstab(df_n_masked['dropout_bin'], df_n_masked['het_bin'])
    # Convert to long format
    ct_long = ct.reset_index().melt(id_vars='dropout_bin',
                                    var_name='het_bin',
                                    value_name='count')
    fig = px.bar(
        ct_long,
        x="dropout_bin",
        y="count",
        color="het_bin",
        color_discrete_sequence=px.colors.sequential.Viridis_r, 
    )
    fig.update_layout(
        xaxis_title=f"Number of dropout masked positions per sample",
        yaxis_title="Number of samples",
        # bargap=0.05,
        legend_title="Heterozygous sites",
        width=1920,
        height=1080,
    )
    fig.update_layout(
        legend=dict(
            x=0.80,  # horizontal position: 0=left, 1=right
            y=0.95,  # vertical position: 0=bottom, 1=top
            xanchor='center',  # align the legend box by its center
            yanchor='top',
            # orientation='h',
            bgcolor='rgba(255,255,255,0.5)',  # semi-transparent background
            bordercolor='black',
            borderwidth=1
        ),
        title=dict(
            text="Distribution of the ~5M samples binned per number of dropout masked positions,<br>stacked per number of heterozygous sites",
            x=0.5,         # 0 = left, 0.5 = center, 1 = right
            xanchor='center',  # Align title's anchor point
            font=dict(size=16),
        ),
    )

    fig_dir = Path(f"{data_dir}/figs/stacked_dropout_masked_het_sites_{param_term}.png")
    fig_dir.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(fig_dir, scale=2)
    
def plot_reference_tree_representativity_stacked_barplot(reference_tree_samples_list, metadata_tsv, data_dir, param_term):
    """Plot the stacked bar plot showing the lineages binned by the number of samples they contain, stacked by two categories: presence or absence in the reference tree.
    """
    with smart_open(metadata_tsv, "rt") as f:
        metadata_df = pd.read_csv(f, sep="\t", low_memory=False)
    metadata_df = metadata_df[['Run', 'Viridian_pangolin_1.29']].copy()
    number_of_unique_lineages = metadata_df['Viridian_pangolin_1.29'].nunique()
    print(f"{number_of_unique_lineages} unique lineages in the initial dataset")
    
    metadata_ref_df = metadata_df.loc[metadata_df['Run'].isin(clean_samples_name_list)].copy()
    n_ref_lineages = metadata_ref_df['Viridian_pangolin_1.29'].nunique()
    print(f"{n_ref_lineages} in the reference tree")
    
    # Try to make a stacked bar plot binning lineages per number of samples that they contain and showing how much of those lineages are in the ref
    # So each bar would be cut in two, showing the number of lineages in this current bin that are in the reference tree (this bin can correspond to the lineages that have 0 - 1000 samples for instance)

    # Count number of samples per lineage
    lineage_counts = metadata_df.groupby('Viridian_pangolin_1.29')['Run'].count().reset_index()
    lineage_counts.columns = ['lineage', 'n_samples']

    # Mark whether the lineage is in the reference tree
    lineage_counts['in_ref'] = lineage_counts['lineage'].isin(metadata_ref_df['Viridian_pangolin_1.29']).map({True: 'In reference', False: 'Not in reference'})

    # Define bins for number of samples per lineage
    bins = [1, 10, 30, 50, 100, 250, 500, 1000, np.inf]
    labels = [f"{bins[i]} - {bins[i+1]-1}" for i in range(len(bins)-2)] + [f"{bins[-2]}+"]
    # labels = [str(bin_int) for bin_int in bins[:-2]] + [f"{bins[-2]}+"]
    lineage_counts['sample_bin'] = pd.cut(lineage_counts['n_samples'], bins=bins, labels=labels, right=False)

    # Count number of lineages per bin and ref status
    ct = lineage_counts.groupby(['sample_bin', 'in_ref'], observed=False).size().reset_index(name='n_lineages')

    # Stacked bar plot
    colours = [px.colors.sequential.Viridis[-6], px.colors.sequential.Viridis[-9]]
    fig = px.bar(
        ct,
        x='sample_bin',
        y='n_lineages',
        color='in_ref',
        text='n_lineages',
        color_discrete_sequence=colours, 
    )

    fig.update_traces(textposition='outside')
    fig.update_layout(
        xaxis_title="Number of samples per lineage",
        yaxis_title="Number of lineages",
        legend_title_text="",
        legend_traceorder="reversed",
        barmode='stack'
    )
    fig.update_layout(
        legend=dict(
            x=0.5,  # horizontal position: 0=left, 1=right
            y=0.93,  # vertical position: 0=bottom, 1=top
            xanchor='center',  # align the legend box by its center
            yanchor='top',
            # orientation='h',
            bgcolor='rgba(255,255,255,0.5)',  # semi-transparent background
            bordercolor='black',
            borderwidth=1,
        ),
        legend_title=dict(
            text="",
        ),
        title=dict(
            text='Distribution of the 3,631 Pango lineages binned per number of samples,<br>stacked by reference tree presence',
            x=0.5,         # 0 = left, 0.5 = center, 1 = right
            xanchor='center',  # Align title's anchor point
            font=dict(size=16),
        ),
    )
    fig_dir = Path(f"{data_dir}/figs/stacked_reference_tree_representativity_{param_term}.png")
    fig_dir.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(fig_dir, scale=2)

if __name__ == "__main__":
    """This script concatenates the various files generated by the first (batched) script
    """
    batches_folder_path = Path(f"{data_dir}/1/")

    # CONCATENATE MUTATIONS FILES
    final_mut_path = Path(f"{data_dir}/2/n_masked_and_masked_mut/masked_mut_{param_term}.tsv")
    # Generate the parent directory if it doesn't exist
    final_mut_path.parent.mkdir(parents=True, exist_ok=True)
    
    if not final_mut_path.exists():
        with open(final_mut_path, "w") as outfile:
            outfile.write("sample_name\tmasked_mutations_masked\tmasked_mutations_random\n")

        with open(final_mut_path, "a") as outfile:
            for mut_tsv_file in batches_folder_path.glob(f"maple_mutations_batch*_{param_term}.tsv"):
                with open(mut_tsv_file, "r") as infile:
                    shutil.copyfileobj(infile, outfile)
                mut_tsv_file.unlink()

    # CONCATENATE N_MASKED FILES
    final_path_n_masked = Path(f"{data_dir}/2/n_masked_and_masked_mut/n_masked_{param_term}.tsv")
    if not final_path_n_masked.exists():
        with open(final_path_n_masked, 'w') as outfile:
            outfile.write("sample_name\tn_masked_cons\tn_masked_masked\tn_het_sites\n")
        with open(final_path_n_masked, 'a') as outfile:
            for n_masked_tsv_file in tqdm(batches_folder_path.glob(f"n_masked_batch_*_{param_term}.tsv"), mininterval=30, desc="Processing n_masked files"):
                with open(n_masked_tsv_file, 'r') as infile:
                    shutil.copyfileobj(infile, outfile)
                n_masked_tsv_file.unlink()

    plot_reference_samples_stacked_barplot(n_masked_path_tsv=final_path_n_masked,
                                           data_dir=data_dir,
                                           param_term=param_term)


    # CONCATENATE ALIGNMENT FILES
    final_clean_tree_path = f"{data_dir}/2/alignment_files/clean_tree_alignment_file_{param_term}.maple"
    maple_alignment_unmasked_file_path = f"{data_dir}/2/alignment_files/unmasked_alignment_{param_term}.maple"
    maple_alignment_random_file_path = f"{data_dir}/2/alignment_files/random_alignment_{param_term}.maple"
    maple_alignment_masked_file_path = f"{data_dir}/2/alignment_files/masked_alignment_{param_term}.maple"
    # Generate parent directory if it doesn't exist
    Path(final_clean_tree_path).parent.mkdir(parents=True, exist_ok=True)
    # For each file, generate it if it does not exist
    if not Path(maple_alignment_unmasked_file_path).exists() and not Path(maple_alignment_unmasked_file_path + '.gz').exists():
        build_maple_file(path_ref_seq, maple_alignment_unmasked_file_path)
    if not Path(maple_alignment_random_file_path).exists() and not Path(maple_alignment_random_file_path + '.gz').exists():
        build_maple_file(path_ref_seq, maple_alignment_random_file_path)
    if not Path(maple_alignment_masked_file_path).exists() and not Path(maple_alignment_masked_file_path + '.gz').exists():
        build_maple_file(path_ref_seq, maple_alignment_masked_file_path)
    if not Path(final_clean_tree_path).exists(): build_maple_file(path_ref_seq, final_clean_tree_path)

    clean_samples_maple_list = []
    current_sample_entry = []
    n_samples = 0
    for clean_tree_file in batches_folder_path.glob(f"clean_tree_batch*_{param_term}.maple"):
        with open(clean_tree_file, "r") as f:
            # Add each entry as a list to the clean_samples_maple_list
            for line in f:
                if line.startswith(">"):
                    n_samples+=1
                    # Save sample_name
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
    unique_samples_list = list(set(clean_samples_maple_list))
    
    if len(unique_samples_list) > 0:
        with open(final_clean_tree_path, "a") as new_f:
            new_f.writelines(unique_samples_list)
        print(f"Number of clean samples after removing duplicates: {len(unique_samples_list)}")
        print(f"Clean tree file saved")
    else:
        print("No clean_tree_files, script is probably being reran")
        
    for cons_align_file in batches_folder_path.glob(f"unmasked_alignment_batch*_{param_term}.maple"):
        # Copy the file content to the MAPLE input folder
        with open(cons_align_file, "r") as f:
            with open(maple_alignment_unmasked_file_path, "a") as new_f:
                shutil.copyfileobj(f, new_f)
        cons_align_file.unlink()

    for random_align_file in batches_folder_path.glob(f"random_alignment_batch*_{param_term}.maple"):
        # Copy the file content to the MAPLE input folder
        with open(random_align_file, "r") as f:
            with open(maple_alignment_random_file_path, "a") as new_f:
                shutil.copyfileobj(f, new_f)
        random_align_file.unlink()

    for masked_align_file in batches_folder_path.glob(f"masked_alignment_batch*_{param_term}.maple"):
        # Copy the file content to the MAPLE input folder
        with open(masked_align_file, "r") as f:
            with open(maple_alignment_masked_file_path, "a") as new_f:
                shutil.copyfileobj(f, new_f)
        masked_align_file.unlink()
    print(f"Alignment files saved")
    
    if compress and not Path(maple_alignment_unmasked_file_path + '.zst').exists():
        compress_file(maple_alignment_unmasked_file_path, threads=num_cores)
    if compress and not Path(maple_alignment_random_file_path + '.zst').exists():
        compress_file(maple_alignment_random_file_path, threads=num_cores)
    if compress and not Path(maple_alignment_masked_file_path + '.zst').exists():
        compress_file(maple_alignment_masked_file_path, threads=num_cores)
    
    # Count number of samples in clean tree file
    clean_samples_name_list = [] 
    n_samples = 0
    with open(final_clean_tree_path, "r") as f:
        next(f)
        next(f)
        for line in f:
            if line.startswith(">"):
                n_samples += 1
                sample_name = line[1:].strip()
                clean_samples_name_list.append(sample_name)
    print(f"Number of samples in clean tree file: {n_samples}")

    plot_reference_tree_representativity_stacked_barplot(reference_tree_samples_list=clean_samples_name_list,
                                                         metadata_tsv=path_metadata_tsv,
                                                         data_dir=data_dir,
                                                         param_term=param_term)


    # GENERATE ALIGNMENT FILES
    # That is, append the reference sequences to the batch alignment files of samples to process through the pipeline
    with open(final_clean_tree_path, "r") as f:
        next(f)
        next(f) # Skip the first two lines (reference sequence)
        clean_tree_content = f.read().encode()

    align_files = list(batches_folder_path.glob(f"maple_alignment_batch*_{param_term}*"))
    print(f"Appending clean samples to {len(align_files)} alignment files using {num_cores} CPU cores...")
    with Pool(processes=num_cores) as pool:
        pool.map(partial(append_clean_tree, content=clean_tree_content), align_files)    
    print("Added clean samples to alignment files")


    # BUILD MUTATION DICTIONARIES
    var_mut_dict = build_mut_dict(final_clean_tree_path) # ref_seq is read from the file
    print(f"Built mutation dictionary with {len(var_mut_dict)} unique variant entries.")
    
    # Build mut_to_variant_dict, where keys are mutations and values are lists of variants having that mutation
    # For simplification and better memory efficiency, we store those as TSV files and recreate them back when needed.
    # The TSV has to have only two columns, unnamed, col 0 are keys, col 1 are values
    # Create dict folder if it doesn't exist
    dict_folder = Path(f"{data_dir}/2/dict/")
    dict_folder.mkdir(parents=True, exist_ok=True)
    dict_path = Path(f"{data_dir}/2/dict/clean_samples_var_to_mut_dict_{param_term}.pickle")
    inv_dict_path = Path(f"{data_dir}/2/dict/clean_samples_mut_to_var_dict_{param_term}.pickle")
    
    if not dict_path.exists() and not inv_dict_path.exists():
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

        # Save variant:mutation dictionary as TSV
        save_pickle_dict(var_mut_dict, dict_path, compress)
        print(f"Variant to mutation dictionary saved to {dict_path}.")

        # Save mutation:variant dictionary as TSV
        save_pickle_dict(mut_to_variant_dict, inv_dict_path, compress)
        print(f"Mutation to variant dictionary saved to {inv_dict_path}.")

    