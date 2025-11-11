# DETECTION OF POTENTIALLY CONTAMINATED AREAS

import platform
from tqdm import tqdm
import gzip
from statistics import median
import lzma
import argparse
import shutil
from pathlib import Path
import random
from hmmlearn import hmm
import numpy as np
from _aux_functions import generate_sh_param_file, compress_file, build_maple_entry, build_maple_file, update_params_file
import signal
from functools import wraps
import time

parser = argparse.ArgumentParser(description='Detect potential contaminated areas depending on several thresholds applied on the Viridian samples.')

parser.add_argument('--n_batch', type=int,
                    help='Number of batches to split the samples into.')
parser.add_argument('--batch_id', type=int,
                    help='Batch id to process.')
parser.add_argument('--data_dir', type=str,
                    help='Directory containing all the files generated during the pipeline.')
parser.add_argument('--samples_dir', type=str,
                    help='Directory containing the sample names and paths to their Viridian output folders.')
parser.add_argument('--param_path', type=str,
                    help='Path to the parameters file.')
parser.add_argument('--depth_thr', type=float,
                    help='Depth threshold for masking.')
parser.add_argument('--max_n_het_sites', type=int,
                    help='Maximum number of heterozygous sites to consider a sample clean enough to be part of the clean tree.')
parser.add_argument('--path_ref_seq', type=str,
                    help='Path to the reference sequence used for the pipeline.')
parser.add_argument('--het_thr', type=float,
                    help='Heterozygosity threshold used to consider a site heterozygous or not.')

args = parser.parse_args()

# Get the arguments
n_batch = args.n_batch
batch_id = args.batch_id
data_dir = args.data_dir
samples_dir = args.samples_dir
param_path = args.param_path
depth_thr = args.depth_thr
max_n_het_sites = args.max_n_het_sites
path_ref_seq = args.path_ref_seq
het_thr = args.het_thr

param_term = f"{depth_thr}_{max_n_het_sites}"

masking_method="thr"


def generate_batchs_pairs(n_batch, path_vdn = samples_dir):
    """Generate a 2D list with structure such that L[i][j] = sample j of batch i.
    Samples are pairs (path, read_name) where path is the path to the qc file.

    Args:
        n_batch (int): Number of batches to generate.
    """
    # Generate a list of ids to chose samples
    batchs_samples_list = [[] for _ in range(n_batch)]

    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        for i, line in enumerate(f):
            if i == 0:
                continue
            # Get the read_name and path
            read_name, path = line.strip().split("\t")
            path = path + "/qc.tsv.gz"
            batchs_samples_list[i % n_batch].append((path, read_name))
            
            
        # for i, batch in tqdm(enumerate(batchs_samples_list), desc="Generating batches", total=n_batch):
        #     with lzma.open(f"{data_storing_dir}/batches/batch_{i}.pkl.xz", "wb") as f:
        #         pickle.dump(batch, f)
        # Not enough memory to store all batches
    
    return batchs_samples_list

def build_masked_pos_dict():
    
    # Build filtering dict
    maskedPoss = {i: False for i in range(1, 29904)}
    # List of 29903 positions with True if the position is masked
    for i in range(1, 73): # Masking the first 72 positions
        maskedPoss[i] = True
        
    for i in range(29769, 29903+1): # Masking the last 134 positions
        maskedPoss[i] = True
        
    maskedPoss[25202]=True # high het
    maskedPoss[21987]=True # high het
    maskedPoss[23118]=True # high het
    maskedPoss[8835]=True # high het
    maskedPoss[15521]=True # high het
    maskedPoss[27507]=True
    maskedPoss[26766]=True
    maskedPoss[8008]=True
    maskedPoss[8012]=True
    maskedPoss[15510]=True
    maskedPoss[17259]=True
    maskedPoss[19413]=True
    maskedPoss[22786]=True
    maskedPoss[22882]=True
    maskedPoss[23948]=True
    maskedPoss[8826]=True
    maskedPoss[8829]=True
    maskedPoss[15854]=True
    maskedPoss[19672]=True
    maskedPoss[21650]=True
    maskedPoss[25296]=True
    maskedPoss[25324]=True
    maskedPoss[25336]=True
    maskedPoss[29687]=True
    maskedPoss[22026]=True
    maskedPoss[22027]=True
    maskedPoss[22028]=True
    maskedPoss[22029]=True
    maskedPoss[22030]=True
    maskedPoss[22031]=True
    maskedPoss[22032]=True
    maskedPoss[22033]=True
    maskedPoss[22034]=True
    maskedPoss[22195]=True
    maskedPoss[22197]=True
    maskedPoss[22198]=True
    maskedPoss[22202]=True
    maskedPoss[22204]=True
    maskedPoss[274]=True
    maskedPoss[4321]=True
    maskedPoss[26530]=True
    maskedPoss[28245]=True
    maskedPoss[28247]=True
    maskedPoss[28249]=True
    maskedPoss[28253]=True
    maskedPoss[28251]=True
    maskedPoss[28254]=True
    
    return maskedPoss


def apply_maple_masking(seq, ref_seq=None):
    # This masking can either be replacing with 'n' or replacing with the reference sequence nucleotide.
    # The second option is more efficient
    lseq = len(seq)
    # Build filtering dict
    maskedPoss = build_masked_pos_dict()
    
    new_seq = []
    for i in range(1, lseq+1):
        if maskedPoss[i]:
            
            if ref_seq:
                ref_seq_nt = ref_seq[i-1]
                new_seq.append(ref_seq_nt.lower())
                
            else:
                new_seq.append("n")
    
        else:
            new_seq.append(seq[i-1])
    
    return new_seq

def retrieve_pos_to_mask_hmm(hmm_model, clean_depth_list, sample_name):
    """Setup HMM, fit the list and output the list of positions to mask.

    Args:
        clean_depth_list (list): list of clean depth values.
    """
    X = np.array(clean_depth_list).reshape(-1, 1)
    # Clip the values at 1800
    X = np.clip(X, 0, 1800)
    print(sample_name, flush=True)
    
    def timed_out_predict(hmm_model, X, sample_name):
        hidden_states = hmm_model.predict(X)
        positions_to_mask = [i for i, state in enumerate(hidden_states) if state == 1 or state == 2]
        return positions_to_mask, True

    # Mask positions where the state is 1 or 2 (contaminated region)
    positions_to_mask, is_processed_correctly = timed_out_predict(hmm_model, X, sample_name)

    return positions_to_mask, is_processed_correctly

def apply_masking(sample_list,
                path_ref_seq,
                path_write_file,
                het_thr,
                depth_thr,
                prop_under_depth_thr, 
                typical_depth_thr,
                n_masked_thr,
                max_n_het_sites,
                param_term,
                data_dir,
                hmm_model=None):
    """Read the Viridian qc file, apply several masking thresholds and
    generate a MAPLE file with the VDN unmasked sequences, our masked sequences, and random masked sequences.

    Args:
        sample_list (list): list of tuples containing the sample name and the path to the file.
        het_thr (float): threshold for the heterozygous proportion.
        depth_thr (float): threshold for the depth.
        prop_under_depth_thr (float): threshold for the maximal proportion of positions with clean depth below the threshold.
        typical_depth_thr (int): median depth threshold.
        path_ref_seq (str): path to the reference sequence.
        path_write_file (str): path to the output file.
        n_masked_thr (int): threshold for the maximal number of masked positions.
        max_n_het_sites (int): maximum number of heterozygous sites to consider a sample clean enough to be part of the clean tree.
        param_term (str): string to identify the parameters used.
        hmm_model (hmm.GaussianHMM): HMM model to use for masking. If None, no HMM masking is applied.
    Returns:
        str: path to the MAPLE alignment file.
    """
    ref = build_maple_file(path_ref_seq, path_write_file)
    n_samples = len(sample_list)
    path_mut_file = f"{data_dir}/1/maple_mutations_batch{batch_id}_{param_term}.tsv"
    with open(path_mut_file, "w") as f:
        f.write("sample_name\tmasked_mutations_masked\tmasked_mutations_random\n")

    clean_tree_files_path = f"{data_dir}/1/clean_tree_batch{batch_id}_{param_term}.maple"

    maple_alignment_unmasked_file_path = f"{data_dir}/1/unmasked_alignment_batch{batch_id}_{param_term}.maple"
    maple_alignment_random_file_path = f"{data_dir}/1/random_alignment_batch{batch_id}_{param_term}.maple"
    maple_alignment_masked_file_path = f"{data_dir}/1/masked_alignment_batch{batch_id}_{param_term}.maple"

    path_store_n_masked = f"{data_dir}/1/n_masked_batch_{batch_id}_{param_term}.tsv"

    # Iterate over the test_list
    for path, read_name in tqdm(sample_list,
                                desc="Processing samples",
                                total=n_samples,
                                unit="sample",
                                mininterval=90):
        clean_depth_list = []
        cons_depth_list = []
        ref_pos_list = []
        unma_seq = []
        
        # unzip and open the file
        with gzip.open(path, "rt") as f:
            # Read the lines and create the needed values to apply the thresholds
            # List of all the clean_depth, list of het_prop, median_depth
            
            # b'Ref_pos\tRef_nt\tCons_pos\tCons_nt\tMasked_cons_nt\tAmplicon\tPrimer\tMask\tTotal_depth\tClean_depth\tCons_depth\tA\ta\tC\tc\tG\tg\tT\tt\tI\ti\tD\td\tX_A\tX_a\tX_C\tX_c\tX_G\tX_g\tX_T\tX_t\tX_I\tX_i\tX_D\tX_d\n'
            
            # Iterate over the lines in the file
            next(f)  # Skip header line
            for line in f:
                # Decode the line and split it into columns
                line = line.strip().split("\t")
                # Extract the values from the columns
                if line[9].isdigit():
                    # Extract the values from the columns
                    ref_pos = int(line[0])
                    clean_depth = int(line[9])
                    cons_depth = int(line[10])
                    # Lowercase the reference and unmasked sequences
                    cons_nt = line[4].lower() # Consensus is Viridian's consensus
                else:
                    ref_pos = int(line[0])
                    clean_depth = 0
                    cons_depth = 0
                    cons_nt = line[4].lower()
                    # Skip the first and last lines where depths col are filled with .
                
                if (len(ref_pos_list) > 0) and (ref_pos == ref_pos_list[-1]):
                    # Insertion
                    unma_seq[-1] += cons_nt
                
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    clean_depth_list.append(clean_depth)
                    cons_depth_list.append(cons_depth)
        
        # Apply the thresholds
        # 1. Check if median depth is above the threshold
        median_depth = median(clean_depth_list)
        
        if median_depth < typical_depth_thr:
            continue
        
        # 2. Check if the proportion of clean depth below the threshold is above the threshold
        depth_int_thr = median_depth * depth_thr
        prop_under_depth = sum(1 for d in clean_depth_list if d < depth_int_thr) / len(clean_depth_list)
        
        if prop_under_depth > prop_under_depth_thr:
            continue
        
        if len(unma_seq) != 29903:
            # MAPLE needs all the sequences to be of the same length
            # With the insertion management, normally this if condition should never be reached
            print(f"Warning: Sequence length is not 29903 for {read_name} with len {len(unma_seq)}. Skipping this sample.")
            continue
        
        
        # 3. MASKING
        # Retrieve the indexes where the het_prop is under the threshold AND the clean depth is under the threshold
        # Build the masked sequence and random masked sequence
        unma_seq = list(unma_seq)
        original_seq = unma_seq.copy()
        # Masking of the very heterozygous or prone to errors positions according to Nicola's papers 
        unma_seq = apply_maple_masking(unma_seq, ref) # WE INPUT THE REFERENCE TO REPLACE WITH THE REF NT INSTEAD OF 'n'
        
        # Decontaminator masking
        masked_seq = unma_seq.copy()
        random_seq = unma_seq.copy()
        masked_pos = []
        shift = random.randint(1000, 29000)
        
        # Count the number of heterozygous sites
        n_het_sites = 0

        for i, (clean_depth, cons_depth) in enumerate(zip(clean_depth_list, cons_depth_list)):
            # Loop to count het sites and potentially mask positions if using threshold-based masking
            het_prop = 1 if clean_depth == 0 else cons_depth / clean_depth

            if het_prop <= 1 - het_thr:
                n_het_sites += 1
            
            if not hmm_model:
                if clean_depth <= max(depth_int_thr, 50):
                    # 50 is the minimum depth for consensus calling in our version
                    if masked_seq[i] != "n":
                        # We only count the positions that are masked by DPCA and NOT MAPLE or VIRIDIAN
                        masked_seq[i] = "n"
                        random_seq[(i + shift) % 29903] = "n"
                        masked_pos.append(i)

        if hmm_model:
            # Setup model
            positions_to_mask, is_processed_correctly = retrieve_pos_to_mask_hmm(hmm_model, clean_depth_list, read_name)
            if is_processed_correctly:
                for pos in positions_to_mask:
                    if masked_seq[pos] != "n":
                        # We only count the positions that are masked by DPCA and NOT MAPLE or VIRIDIAN
                        masked_seq[pos] = "n"
                        random_seq[(pos + shift) % 29903] = "n"
                        masked_pos.append(pos)
            else:
                # Sample could not be processed correctly by the HMM, we skip it
                print(f"Skipping sample {read_name} as it could not be processed correctly by the HMM.")
                continue
                
        
        # KEEPING SAMPLES TO BUILD THE CLEAN TREE
        if n_het_sites <= max_n_het_sites:
            # print(f"Warning: No masked positions for {read_name}. Storing this sample for the clean tree.")
            
            with open(clean_tree_files_path, "a") as f:
                # SHOULD I USE unma_seq OR ORIGINAL SEQ? unma_seq has the MAPLE masking while original_seq is the exact Viridian consensus
                cons_entry = build_maple_entry(unma_seq, ref, f"{read_name}")
                f.write(f"{cons_entry}\n")
                
            # This sample cannot be tested for contamination as it is used to build the clean tree
            continue
        
        # 5. Check the number of masked positions and apply the threshold
        n_masked_masked = sum(1 for c in masked_seq if c == "n")
        n_random_masked = sum(1 for c in random_seq if c == "n")
        n_masked_cons = sum(1 for c in unma_seq if c == "n")
        n_het_sites = int(n_het_sites)
        
        # Store these values in a tsv file
        with open(path_store_n_masked, "a") as f:
            f.write(f"{read_name}\t{n_masked_cons}\t{n_masked_masked}\t{n_het_sites}\n")
        
        if n_masked_masked > n_masked_thr:
            continue
        
        n_diff = n_masked_masked - n_random_masked
        # We want n_diff to be 0 after this loop
        while n_diff != 0:
            # print(n_diff, n_masked_masked, n_random_masked)
            # Mask n_diff consecutive positions in the random sequence
            start = random.randint(100, 29800)
            for i in range(n_diff):
                random_seq[(start + i) % 29903] = "n"
            
            n_random_masked = sum(1 for c in random_seq if c == "n")
            n_diff = n_masked_masked - n_random_masked
        
        # 6. Add the VDN unmasked consensus, our masked sequence and the random masked sequence to the output file
        with open(path_write_file, "a") as f:
            cons_entry = build_maple_entry(unma_seq, ref, f"unmasked_{read_name}")
            masked_entry = build_maple_entry(masked_seq, ref, f"masked_{read_name}")
            random_entry = build_maple_entry(random_seq, ref, f"random_{read_name}")
            f.write(f"{cons_entry}\n{masked_entry}\n{random_entry}\n")
            
        with open(maple_alignment_unmasked_file_path, "a") as f:
            f.write(f"{cons_entry}\n")
        with open(maple_alignment_random_file_path, "a") as f:
            f.write(f"{random_entry}\n")
        with open(maple_alignment_masked_file_path, "a") as f:
            f.write(f"{masked_entry}\n")
        
        
        # 7. Store the mutations masked in the masked and random sequences
        cons_mut = []
        masked_seq_mut = []
        random_seq_mut = []
        for line in cons_entry.split("\n")[1:]:
            if line.startswith("n") or line.startswith("-"):
                continue
            parts = line.split("\t")
            nuc = parts[0]
            pos = int(parts[1])
            ref_nt = ref[pos - 1]
            # Higher case for the nt
            nuc = nuc.upper()
            ref_nt = ref_nt.upper()
            cons_mut.append(f"{ref_nt}{pos}{nuc}")
            
        for line in masked_entry.split("\n")[1:]:
            if line.startswith("n") or line.startswith("-"):
                continue
            parts = line.split("\t")
            nuc = parts[0]
            pos = int(parts[1])
            ref_nt = ref[pos - 1]
            nuc = nuc.upper()
            ref_nt = ref_nt.upper()
            masked_seq_mut.append(f"{ref_nt}{pos}{nuc}")
            
        for line in random_entry.split("\n")[1:]:
            if line.startswith("n") or line.startswith("-"):
                continue
            parts = line.split("\t")
            nuc = parts[0]
            pos = int(parts[1])
            ref_nt = ref[pos - 1]
            nuc = nuc.upper()
            ref_nt = ref_nt.upper()
            random_seq_mut.append(f"{ref_nt}{pos}{nuc}")
            
        masked_mut = set(cons_mut) - set(masked_seq_mut)
        rand_masked_mut = set(cons_mut) - set(random_seq_mut)
     
        # Write the mutations to a file
        with open(path_mut_file, "a") as f:
            masked_mut = ';'.join(masked_mut)
            masked_mut = masked_mut.upper()
            
            rand_masked_mut = ';'.join(rand_masked_mut)
            rand_masked_mut = rand_masked_mut.upper()
            f.write(f"{read_name}\t{masked_mut}\t{rand_masked_mut}\n")

    return path_write_file


if __name__ == "__main__":
    # Generate the batches
    print(platform.python_implementation(), platform.python_version())
    print(f"Data directory: {data_dir}")
    # Generate the directory to store the batches if it does not exist
    
    if batch_id==1:
        
        # Update the params.py file with data_dir
        updates = {
            "n_batch": n_batch,
            "data_dir": data_dir,
            "samples_dir": samples_dir,
            "param_path": param_path,
            "depth_thr": depth_thr,
            "max_n_het_sites": max_n_het_sites,
            "path_ref_seq": path_ref_seq,
            "het_thr": het_thr
        }

        update_params_file(param_path, updates)
        print("Generated params.py file.")
        param_sh_path = param_path.replace(".py", ".sh")
        generate_sh_param_file(param_sh_path)
        # Remove the files in the batches folder if they exist
        batches_folder = Path(f"{data_dir}/1/")
        batches_folder.mkdir(parents=True, exist_ok=True)
        for file in batches_folder.glob(f"*{param_term}*"):
            file.unlink()
        print("Cleared the batches folder.")

    else:
        print("Waiting for process 1 to clear the batches folder...")
        time.sleep(90)

    batchs_samples_list = generate_batchs_pairs(n_batch)
    
    sample_list = batchs_samples_list[batch_id]
    
    del batchs_samples_list
    
    print(param_term)
    path_write_file = f"{data_dir}/1/maple_alignment_batch{batch_id}_{param_term}.maple"

    if masking_method == "hmm":
        # Initialize the HMM model with the parameters from params.py
        hmm_model = hmm.GaussianHMM(
        n_components=3,
        covariance_type="diag",
        min_covar=1e-3,
        startprob_prior=1.0,
        transmat_prior=1.0,
        means_prior=0,
        means_weight=0,
        covars_prior=1e-2,
        covars_weight=1,
        algorithm="viterbi",
        random_state=7,
        n_iter=20,
        tol=1e-2,
        verbose=False,
        )
        
        hmm_model.means_ = np.array(means).reshape(-1, 1)
        hmm_model.covars_ = np.array(covars).reshape(-1, 1)
        hmm_model.startprob_ = np.array(startprob, dtype=float).ravel()
        hmm_model.transmat_  = np.array(transmat, dtype=float).reshape(3, 3)

    else:
        hmm_model = None
        
    prop_under_depth_thr = 1
    # unused
    typical_depth_thr = 0
    # unused
    n_masked_thr = 30000
    # unused

    # Run the function on each batch
    path_write_file = apply_masking(sample_list,
                                    het_thr=het_thr,
                                    depth_thr=depth_thr,
                                    prop_under_depth_thr=prop_under_depth_thr,
                                    typical_depth_thr=typical_depth_thr,
                                    n_masked_thr=n_masked_thr,
                                    max_n_het_sites=max_n_het_sites,
                                    path_ref_seq=path_ref_seq,
                                    path_write_file=path_write_file,
                                    param_term=param_term,
                                    hmm_model=hmm_model,
                                    data_dir=data_dir)

    # Compress the file
    # compress_file(path_write_file)
    
    # Save the result
    print(f"Batch {batch_id} done")
    
    # Save a done file
    done_file_path = Path(f"{data_dir}/done_files/1_maple_alignment_batch_{batch_id}.done")
    done_file_path.parent.mkdir(parents=True, exist_ok=True)
    done_file_path.touch()
    
