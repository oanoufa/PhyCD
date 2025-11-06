# This script should go over the tsv results of MAPLE sample placement and
# generate the distributions of the branch lengths for the three types of sequences
# (viridian, masked, random).

# The goal is to see whether the masked samples are closer to the tree than the viridian ones or not.

import pickle
from pathlib import Path
import gzip
import lzma
from tqdm import tqdm
import pandas as pd
from collections import Counter
import _params
from _aux_functions import update_params_file
import argparse
import csv


parser = argparse.ArgumentParser(description='Process the concatenated output of the MAPLE sample placement.')

parser.add_argument('--masked_or_random', type=str,
                    help='Either masked or random, to process both in parallel')

args = parser.parse_args()

# Get the arguments
masked_or_random = args.masked_or_random

path_ref_seq = _params.path_ref_seq
inputTree = _params.inputTree
inputRates = _params.inputRates
n_diff_mut = _params.n_diff_mut
masked_max_dist = _params.masked_max_dist
masking_ratio = _params.masking_ratio
param_term = _params.param_term


def remove_samples_without_three_types(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove samples that do not have all three types of sequences (unmasked, masked, random).
    """
    # Count the number of unique sample types for each read name
    counts = df.groupby("sample_name")["sample_type"].nunique()
    
    # Filter to keep only read names with all three types
    invalid_read_names = counts[counts < 3].index
    
    # Filter the original DataFrame to keep only those read names
    df = df[~df["sample_name"].isin(invalid_read_names)]
    
    df = df.sort_values(by=["sample_name"])
    df = df.reset_index(drop=True)
    
    return df

def remove_samples_placed_on_self(row):
    
    sample_name = row['sample_name']
    cons_placement = row['unmasked_placement']
    cons_placement = cons_placement.split(":")[0]  # Get the placement without support
    
    if sample_name == cons_placement:
        # Remove the row
        return True
    return False

def generate_sample_list(sample_names):
    """Get the sample names in the result dataset and retrieve their list of heterozygosity values.
    """
    path_vdn = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"
    # Generate a list of ids to chose samples
    batchs_samples_list = []

    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        next(f)  # Skip header line
        for line in tqdm(f, desc="Generating sample list",
                            mininterval=120):
            # Get the read_name and path
            read_name, path = line.strip().split("\t")
            path = path + "/qc.tsv.gz"
            if read_name in sample_names:
                batchs_samples_list.append((path, read_name))
                sample_names.remove(read_name)
            if len(sample_names) == 0:
                break
    return batchs_samples_list

def build_het_sites_dict(sample_list, het_thr):
    """Take as input a set of sample names, find their qc file,
    iterate over it to store all the heterozygous sites, the nt of the minor allele and its proportion.
    The final dictionary should look like this:
    {
        sample_name1: []
            (pos1, minor_nt, prop),
            (pos2, minor_nt, prop),
            ...
        ],
        sample_name2: [
            (pos1, minor_nt, prop),
            (pos2, minor_nt, prop),
            ...
        ],
        ...
    }

    Args:
        sample_list (set): Set of sample names to process
    """
    # Initialize the het_sites_dict with the keys being the sample names
    het_sites_dict = {sample_name: [] for sample_name in sample_list}
    print(f"Generating sample list: {len(sample_list)} {masked_or_random} samples to process")
    name_path_list = generate_sample_list(sample_list)
    print(f"Sample list generated: {len(name_path_list)} samples to process")

    for path, read_name in tqdm(name_path_list,
                                desc="Processing samples",
                                unit="sample",
                                mininterval=120):
        
        ref_pos_list = []
        unma_seq = []
        
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
                    
                    # find minor allele (second highest)
                    max_idx = counts.index(max(counts))
                    counts[max_idx] = 0
                    minor_idx = counts.index(max(counts))
                    minor_allele = "ACGT"[minor_idx]
                    het_depth = counts[minor_idx] / clean_depth
                    
                    if het_depth >= het_thr and clean_depth >= 50:
                        het_sites_dict[read_name].append((ref_pos, minor_allele, het_depth))
                    
                    cons_nt = line[4].lower()
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

    return het_sites_dict

if __name__ == "__main__":
    
    print(f"Processing {masked_or_random} samples...")

    # Define the path to the folder containing the tsv files
    path_store = Path(f"/nfs/research/goldman/anoufa/data/MAPLE_output/processed_placements/processed_placements_results_{masked_or_random}_{param_term}.tsv")
    
    # If path_store already exists, add a suffix to avoid overwriting
    if path_store.exists():
        suffix = 1
        while path_store.with_suffix(f".{suffix}.tsv").exists():
            suffix += 1
        path_store = path_store.with_suffix(f".{suffix}.tsv")

    with open(path_store, "w") as f:
        f.write("sample_type\tsample_name\tsampleBlength\n")

    unit = 1/29903
    mutations_dict = {}
    unmasked_mut_dict = {}
    masked_mut_dict = {}
    unmasked_placement_dict = {}
    masked_placement_dict = {}

    sample_placement_output_path = Path(f"/nfs/research/goldman/anoufa/data/MAPLE_output/output_FULL_metaData_samplePlacements_{param_term}.tsv.gz")
    # Read the tsv file
    
    with open(path_store, "a") as f_out:

        with gzip.open(sample_placement_output_path, "rt") as file:

            reader = csv.reader(file, delimiter="\t")
            header = next(reader)  # skip header

            # Extract the sample lengths for the three types of sequences
            for line in tqdm(reader, desc=f"Processing {sample_placement_output_path.name}", unit="lines", mininterval=30):
                # Columns:  index   sample	                placements	                    optimizedBlengths (<topBlength>/<bottomBlength>/<sampleBlength>)   mutations
                # Example line: 0	unmasked_ERR4239172	in2495143:0.9997103398863393	in2495143:(6.68997107947098e-05/0/0)	T10029C;G21618C;G22917T;A22995C;T23063A;A23604...
                
                # Skip empty lines
                if line[1] == '':
                    continue
                # We are interested in the sampleBlengths of the placement with highest support.
                sample = line[0]
                placements = line[1]
                optimizedBlengths = line[2]
                mutations = line[3]
                
                sample_type, sample_name = sample.split("_", 1)  # Split the sample name to get the type and read name
                
                if "empty possiblePlacements" in placements:
                    continue
                placements = placements.split(";")[0]  # Get the first placement (highest support)


                optimizedBlengths = line[2]
                sampleBlength = optimizedBlengths.split(";")[0] # Get the first placement (highest support)
                sampleBlength = sampleBlength.split("/")[2] # We get this part: 0)
                sampleBlength = sampleBlength[:-1] # Remove the last character (the closing parenthesis)
                sampleBlength = float(sampleBlength)
                unit_sampleBlength = sampleBlength / unit

                
                # Handle O(0.000000/0.000001/0.000772/0.999227)22942A mutations
                mutations = mutations.split(";")
                for i, mut in enumerate(mutations):
                    if "O(" in mut:
                        probas = mut.split("/")
                        A_prob = float(probas[0][2:])
                        C_prob = float(probas[1])
                        G_prob = float(probas[2])
                        T_prob = float(probas[3].split(")")[0])
                        
                        prob_dict = { "A": A_prob, "C": C_prob, "G": G_prob, "T": T_prob }
                        rest = probas[3].split(")")[1]
                        dir_nt = rest[-1]
                        dir_prob = prob_dict[dir_nt]
                        if dir_prob > 0.001:
                            mut = None
                        else:
                            max_nt = max(prob_dict, key=prob_dict.get)
                            mut = max_nt + rest
                        
                        mutations[i] = mut
                mutations = ";".join([m for m in mutations if m is not None])
                if mutations == "":
                    mutations = None
                
                # Retrieve the placements for the unmasked and masked samples
                if sample_type == "unmasked":
                    unmasked_placement_dict[sample_name] = placements
                    unmasked_mut_dict[sample_name] = mutations
                
                elif sample_type == masked_or_random:
                    masked_placement_dict[sample_name] = placements
                    masked_mut_dict[sample_name] = mutations
                
                new_line = f"{sample_type}\t{sample_name}\t{unit_sampleBlength}\n"
                f_out.write(new_line)                                    
                    
    print('Pandas filtering...', flush=True)
    df = pd.read_csv(path_store, sep="\t")
    print(f"Initial number of rows: {len(df)}", flush=True)
    df = remove_samples_without_three_types(df)
    print(f"Number of rows after removing samples without three types: {len(df)}", flush=True)
    
    # Pivot to get one row per sample_name with columns for each sample type
    df = df.pivot(index='sample_name', columns='sample_type', values='sampleBlength')
    df.reset_index(inplace=True)
    df.columns.name = None  # Remove the columns name
    print(f"Number of rows after pivoting the df: {len(df)}", flush=True)
    print(f"Sample types in the df: {df.columns.tolist()}", flush=True)
    
    print(f"Median and mean branch lengths for unmasked samples: {df['unmasked'].median()} and {df['unmasked'].mean()}", flush=True)

    # Filtering interesting samples:
    # Remove samples where unmasked_dist is under 1 (SAMPLES ALREADY WELL PLACED)
    df = df[df['unmasked'] >= 1]
    print(f"Number of rows after removing samples with unmasked dist under 1: {len(df)}", flush=True)
    
    # Remove samples where masked_dist is over masked_max_dist
    df = df[df[masked_or_random] <= masked_max_dist]
    print(f"Number of rows after removing samples not respecting masked_max_dist {masked_max_dist}: {len(df)}", flush=True)
    
    # n_diff_mut: Difference in distance to the tree between unmasked and masked samples
    
    # Keep only if (df['unmasked'] >= df['masked'] + n_diff_mut)
    df = df[(df['unmasked'] >= df[masked_or_random] + n_diff_mut)]
    print(f"Number of rows after removing samples not respecting min n_diff_mut {n_diff_mut}: {len(df)}", flush=True)
    
    # min_ratio_prop: Minimum ratio between the proportion of genome masked and the proportion between the unmasked and masked branch lengths:
    # Example: Branch lengths 3 1 --> masked 66% of the problematic mutations, prop of genome masked 30% 
    # If min_ratio_prop == 2 we retain this sample but if min_ratio_prop == 3 we don't keep it
    n_masked_path = Path(f"/nfs/research/goldman/anoufa/data/MAPLE_output/n_masked_and_masked_mut/n_masked_{param_term}.tsv")
    n_masked_df = pd.read_csv(n_masked_path, sep="\t", names=["sample_name", "n_masked_cons", "n_masked_masked", "n_het_sites"], header=None)
    
    print(f"Median and mean positions masked in unmasked samples: {n_masked_df['n_masked_cons'].median()} and {n_masked_df['n_masked_cons'].mean()}", flush=True)
    print(f"Median and mean positions masked in random/masked samples: {n_masked_df['n_masked_masked'].median()} and {n_masked_df['n_masked_masked'].mean()}", flush=True)

    df = pd.merge(df, n_masked_df, on="sample_name", how="left")
    print(f"Number of rows after merging with n_masked_df: {len(df)}", flush=True)

    df['prop_gen_masked'] = (df['n_masked_masked'] - df['n_masked_cons']) / 29903
    df['prop_dist_reduced'] = (df['unmasked'] - df[masked_or_random]) / df['unmasked']
    # df = df[df['prop_dist_reduced'] >= (masking_ratio * df['prop_gen_masked'])]
    print(f"Number of rows after removing samples not respecting masking_ratio {masking_ratio}: {len(df)}", flush=True)

    df.sort_values(by="unmasked", ascending=False, inplace=True)
    
    sample_names = df['sample_name'].tolist()
    
    # Join the df with the mutations and placements dictionaries
    unmasked_mut_df = pd.DataFrame(list(unmasked_mut_dict.items()), columns=['sample_name', 'unmasked_mutations_to_tree'])
    unmasked_mut_df = unmasked_mut_df[unmasked_mut_df['sample_name'].isin(sample_names)]
    masked_mut_df = pd.DataFrame(list(masked_mut_dict.items()), columns=['sample_name', f'{masked_or_random}_mutations_to_tree'])
    masked_mut_df = masked_mut_df[masked_mut_df['sample_name'].isin(sample_names)]
    
    unmasked_placements_df = pd.DataFrame(list(unmasked_placement_dict.items()), columns=['sample_name', 'unmasked_placement'])
    unmasked_placements_df = unmasked_placements_df[unmasked_placements_df['sample_name'].isin(sample_names)]
    # split unmasked_placement in two columns: unmasked_placement and unmasked_support
    unmasked_placements_df[['unmasked_placement', 'unmasked_support']] = unmasked_placements_df['unmasked_placement'].str.split(":", expand=True)
    unmasked_placements_df['unmasked_support'] = unmasked_placements_df['unmasked_support'].astype(float)
    
    masked_placements_df = pd.DataFrame(list(masked_placement_dict.items()), columns=['sample_name', f'{masked_or_random}_placement'])
    masked_placements_df = masked_placements_df[masked_placements_df['sample_name'].isin(sample_names)]
    # split masked_placement in two columns: masked_placement and masked_support
    masked_placements_df[[f'{masked_or_random}_placement', f'{masked_or_random}_support']] = masked_placements_df[f'{masked_or_random}_placement'].str.split(":", expand=True)
    masked_placements_df[f'{masked_or_random}_support'] = masked_placements_df[f'{masked_or_random}_support'].astype(float)
    
    # Join with df containing the mutations masked during gen_maple_file
    path_df_mut_masked = f"/nfs/research/goldman/anoufa/data/MAPLE_output/n_masked_and_masked_mut/masked_mut_{param_term}.tsv"
    df_mut_masked = pd.read_csv(path_df_mut_masked, sep='\t')
    df_mut_masked = df_mut_masked[df_mut_masked['sample_name'].isin(sample_names)]
    df_mut_masked = df_mut_masked[['sample_name', f'masked_mutations_{masked_or_random}']]
    
    # Merge the dataframes
    print("Merging dataframes...")
    print(f"unmasked mutations df: {unmasked_mut_df.shape}, Masked mutations df: {masked_mut_df.shape}, unmasked placements df: {unmasked_placements_df.shape}, Masked placements df: {masked_placements_df.shape}")
    df = pd.merge(df, unmasked_mut_df, on="sample_name", how="left")
    df = pd.merge(df, masked_mut_df, on="sample_name", how="left")
    df = pd.merge(df, unmasked_placements_df, on="sample_name", how="left")
    df = pd.merge(df, masked_placements_df, on="sample_name", how="left")
    df = pd.merge(df, df_mut_masked, on="sample_name", how="left")
    print(f"FULL df: {df.shape}")
    
    # Remove samples where the highest placement is itself
    df = df[~df.apply(remove_samples_placed_on_self, axis=1)]
    print(f"Number of rows after removing samples placed on themselves: {len(df)}", flush=True)
    
    
    df['type'] = masked_or_random.capitalize()

    df['dist_diff'] = df['unmasked'] - df[masked_or_random]

    df['masking_ratio'] = df['prop_dist_reduced'] / df['prop_gen_masked']

    masked_or_random_inv = 'random' if masked_or_random == 'masked' else 'masked'
    
    df.drop(columns=[masked_or_random_inv], inplace=True)

    df.rename(columns={
        masked_or_random: 'distance',
        'n_masked_masked': 'n_masked',
        'masked_mutations_to_tree': 'mutations_to_tree',
        f'{masked_or_random}_placement': 'placement',
        f'{masked_or_random}_support': 'support',
        f'masked_mutations_{masked_or_random}': 'mutations_masked'}, inplace=True)
    
    
    # Build heterozygous sites dict
    print("Building heterozygous sites dictionary...")
    sample_names = set(df["sample_name"].tolist())
    het_thr = 0.10
    het_sites_dict = build_het_sites_dict(sample_names, het_thr)
    het_sites_path = Path(f"/nfs/research/goldman/anoufa/data/MAPLE_input/dict/het_sites_dict_{masked_or_random}_{param_term}.pickle")
    with open(het_sites_path, "wb") as f:
        pickle.dump(het_sites_dict, f)
    print(f"Heterozygous sites dictionary saved to {het_sites_path}")
    
    # Save df
    df.to_csv(path_store, sep="\t", index=False)
    print(f"Processed placements results saved to {path_store}")
    
    # Add processed_placement to _params.py file and het_sites_dict path
    params_path = "/nfs/research/goldman/anoufa/src/dpca/_params.py"
    update_params_file(params_path,
                       {f'processed_placement_{masked_or_random}':str(path_store), 
                        f'het_sites_dict_{masked_or_random}':str(het_sites_path)}
                       )

    print(f"Processed placements results saved to {path_store}")
        