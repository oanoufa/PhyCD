import pandas as pd
from tqdm import tqdm

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


def look_for_closest_variant(row, variant_mut_dict):
    """
    Optimized version to look for the closest variant.
    Assumes preprocessed_variant_mut_dict has values as sets.
    """
    sample_mut_list = [] # Initialize as empty list
    
    if pd.isna(row['mutations']) or row['mutations'] == "":
        sample_mut_list = row['consensus_mutations'].split(";")
        sample_mut_set = set(sample_mut_list)  # Convert to set for faster intersection later
    else:
        masked_mut_l = row['mutations'].split(";")
        diff_mut_l = row['consensus_mutations'].split(";")
        
        # Take intersection of the two lists using sets
        # This handles both `cons_placement == masked_placement` and `else` cases
        sample_mut_set = set(masked_mut_l) # & set(diff_mut_l)
        # Should I really use the intersection? Or just the masked_mut_l?

    closest_variant = None
    max_closeness = -1
    max_shared_mutations_str = "No shared mutations" # Initialize here

    # --- Step 2: Optimized loop for finding closest variant ---
    # We iterate over the preprocessed dict where values are already sets
    for variant, variant_mut_set in variant_mut_dict.items():
        # Calculate intersection directly with sets
        shared_mutations_set = sample_mut_set & variant_mut_set
        
        mut_in_common = len(shared_mutations_set)
        
        # Avoid division by zero if variant_mut_set is empty (shouldn't happen with real data, but good practice)
        if not variant_mut_set:
            closeness = 0
        else:
            # Goal of this metric is to always take the sample with the most mutations in common and penalize for the number of mutations in the variant
            closeness = mut_in_common - (len(variant_mut_set) / 29903)
        
        # Optimize shared_mutations string construction
        current_shared_mutations_list = sorted(list(shared_mutations_set), key=sorting_key) # Sort here once
        current_shared_mutations_str = ";".join(current_shared_mutations_list)

        # Logic for determining closest_variant and max_closeness
        if closeness > max_closeness and closeness > 0:
            max_closeness = closeness
            closest_variant = variant[1:] # Assuming variant is like ">variant_name"
            max_shared_mutations_str = current_shared_mutations_str
        elif closeness == max_closeness and closeness > 0: # Only concatenate if closeness is positive and ties max_closeness
            # Append new variant name
            if closest_variant is None: # Should not happen if first clause set it
                closest_variant = variant[1:]
            else:
                closest_variant += ", " + variant[1:]
            
            # For simplicity, if they tie and the shared mutations are the same, don't re-add
            if current_shared_mutations_str not in max_shared_mutations_str.split(';'): # Basic check
                 max_shared_mutations_str += ", " + current_shared_mutations_str


    if closest_variant is None: # This catch-all should now be robust
        closest_variant = "No close variant found"
        # max_shared_mutations_str is already "No shared mutations" from initialization

    row['closest_variant'] = closest_variant
    row['closeness'] = max_closeness # Ensure closeness is also set, it was missing in your original snippet
    row['mut_in_common'] = max_shared_mutations_str

    return row

if __name__ == "__main__":
    
    tqdm.pandas(desc="Processing placements", mininterval=60)
    
    placement_result_path = "/nfs/research/goldman/anoufa/data/MAPLE_output/processed_placements_results_4_5.tsv"
    var_mut_d_path = "/nfs/research/goldman/anoufa/data/MAPLE_input/2M_mut_dict_MPL_REF.pickle"
    placement_df = pd.read_csv(placement_result_path, sep="\t")
    variant_mut_dict = pd.read_pickle(var_mut_d_path)
    placement_df = placement_df.progress_apply(look_for_closest_variant, axis=1, variant_mut_dict=variant_mut_dict)
    
    df_aug_path = placement_result_path.replace(".tsv", "_aug.tsv")
    placement_df.to_csv(df_aug_path, sep="\t", index=False)
    print(f"Augmented placement results saved to {df_aug_path}")