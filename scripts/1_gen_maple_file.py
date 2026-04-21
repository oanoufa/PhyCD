# DETECTION OF POTENTIALLY CONTAMINATED AREAS

import platform
from tqdm import tqdm
from statistics import median
import lzma
import argparse
from pathlib import Path
import random
from _aux_functions import build_maple_entry, build_maple_file, compress_file, smart_open
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
parser.add_argument('--depth_thr', type=float,
                    help='Depth threshold for masking.')
parser.add_argument('--max_n_het_sites', type=int,
                    help='Maximum number of heterozygous sites to consider a sample clean enough to be part of the clean tree.')
parser.add_argument('--het_thr', type=float,
                    help='Heterozygosity threshold used to consider a site heterozygous or not.')
parser.add_argument('--max_dropout_masked', type=int,
                    help='Maximum number of dropout masked positions to consider a sample clean enough to be part of the clean tree.')
parser.add_argument('--param_term', type=str,
                    help='str built from the initial parameters used to recognize the files.')
parser.add_argument('--path_ref_seq', type=str,
                    help='Path to the reference sequence used for the pipeline.')
parser.add_argument('--compress', type=int,
                    help='1 or 0, enables compression.')
args = parser.parse_args()

# Get the arguments
n_batch = args.n_batch
batch_id = args.batch_id
data_dir = args.data_dir
samples_dir = args.samples_dir
depth_thr = args.depth_thr
max_n_het_sites = args.max_n_het_sites
het_thr = args.het_thr
max_dropout_masked = args.max_dropout_masked
param_term = args.param_term
path_ref_seq = args.path_ref_seq
compress = args.compress

if compress == 1:
    compress = True
else:
    compress = False

def generate_batchs_pairs(n_batch, path_vdn = samples_dir):
    """Generate a 2D list with structure such that L[i][j] = sample j of batch i.
    Samples are pairs (path, sample_name) where path is the path to the qc file.

    Args:
        n_batch (int): Number of batches to generate.
        path_vdn (str): Path to the file containing the list of (sample_name, path_to_vdn_folder)
    """
    # Generate a list of ids to chose samples
    batchs_samples_list = [[] for _ in range(n_batch)]

    with smart_open(path_vdn, "rt") as f:
        # Each line is composed of sample_name and path to qc file
        # Iterate over the lines and get the path, sample_name if the line's index is in chosen_samples
        for i, line in enumerate(f):
            if i == 0:
                continue
            # Get the sample_name and path
            sample_name, path = line.strip().split("\t")
            path = path + "/qc.tsv.gz"
            batchs_samples_list[i % n_batch].append((path, sample_name))
    
    return batchs_samples_list

def build_masked_pos_dict():
    """Build a dictionary containing all the positions as keys [1, 29903] and True if the position should be masked as value.
    The dictionary is 1-based.

    Returns:
        dict: masked_pos_dict
    """
    
    # Build filtering dict
    maskedPoss = {i: False for i in range(1, 29904)}
    # List of 29903 positions with True if the position is masked
    for i in range(1, 73): # Masking the first 72 positions
        maskedPoss[i] = True
        
    for i in range(29769, 29903+1): # Masking the last 134 positions
        maskedPoss[i] = True

    maskedPoss[274]=True
    maskedPoss[4321]=True
    maskedPoss[8008]=True
    maskedPoss[8012]=True
    maskedPoss[8826]=True
    maskedPoss[8829]=True
    maskedPoss[8835]=True # high het
    maskedPoss[15510]=True
    maskedPoss[15521]=True # high het
    maskedPoss[15854]=True
    maskedPoss[17259]=True
    maskedPoss[19413]=True
    maskedPoss[19672]=True
    maskedPoss[21650]=True
    maskedPoss[21987]=True # high het
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
    maskedPoss[22786]=True
    maskedPoss[22882]=True
    maskedPoss[23118]=True # high het
    maskedPoss[23948]=True
    maskedPoss[25202]=True # high het
    maskedPoss[25296]=True
    maskedPoss[25324]=True
    maskedPoss[25336]=True
    maskedPoss[26530]=True
    maskedPoss[26766]=True
    maskedPoss[27507]=True
    maskedPoss[28245]=True
    maskedPoss[28247]=True
    maskedPoss[28249]=True
    maskedPoss[28251]=True
    maskedPoss[28253]=True
    maskedPoss[28254]=True
    maskedPoss[29687]=True

    # The positions masked are: 
    # 1-72, 29769-29903,
    # 274, 4321, 8008, 8012, 8826, 8829, 8835, 15510, 15521, 15854, 17259, 19413, 19672, 21650, 21987,
    # 22026-22034, 22195, 22197, 22198, 22202, 22204, 22786, 22882, 23118, 23948, 25202, 25296, 25324,
    # 25336, 26530, 26766, 27507, 28245-28247-28249-28251-28253-28254 and 29687.
    
    return maskedPoss

def apply_maple_masking(seq, ref_seq=None):
    """Apply the initial masking as described in the paper.
    We mask the first 72 positions, the last 134 and 47 positions affected by recurrent errors (presented in the function above).
    We either mask by replacing with N or by replacing with the ref nt if ref is given.

    Args:
        seq (list): Sequence as a list
        ref_seq (list, optional): If given, mask by replacing with the reference. Defaults to None.

    Returns:
        list: Sequence masked.
    """
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

def reset_file(path):
    """Truncate the file so it becomes empty."""
    with open(path, "w"):
        pass


def apply_masking(sample_list,
                  path_ref_seq,
                  het_thr,
                  depth_thr,
                  max_dropout_masked,
                  max_n_het_sites,
                  param_term,
                  data_dir,
                  compress,
                  ):
    """Read the Viridian qc file, apply several masking and
    generate a MAPLE file with the dropout unmasked sequences, dropout masked sequences, and randomly masked sequences.
    This function also generates auxiliary data on the samples used in the rest of the pipeline.

    Args:
        sample_list (list): List of (path_to_vdn_qc_file, sample_run_name)
        path_ref_seq (str): Path to the reference sequence used for the analysis
        het_thr (float): Minimal proportion of minor allele to consider a site as heterozygous
        depth_thr (float): Proportion of the median depth under which sites are dropout masked
        max_dropout_masked (int): Max number of dropout masked sites in a filtered sample
        max_n_het_sites (int): Max number of het sites in a filtered sample
        param_term (str): Identifier of the current pipeline
        data_dir (str): Path to the data directory of the current pipeline
    """
    path_write_file = f"{data_dir}/1/maple_alignment_batch{batch_id}_{param_term}.maple"
    ref = build_maple_file(path_ref_seq, path_write_file)
    n_samples = len(sample_list)
    path_mut_file = f"{data_dir}/1/maple_mutations_batch{batch_id}_{param_term}.tsv"
    reset_file(path_mut_file)
    clean_tree_files_path = f"{data_dir}/1/clean_tree_batch{batch_id}_{param_term}.maple"
    reset_file(clean_tree_files_path)
    maple_alignment_unmasked_file_path = f"{data_dir}/1/unmasked_alignment_batch{batch_id}_{param_term}.maple"
    reset_file(maple_alignment_unmasked_file_path)
    maple_alignment_random_file_path = f"{data_dir}/1/random_alignment_batch{batch_id}_{param_term}.maple"
    reset_file(maple_alignment_random_file_path)
    maple_alignment_masked_file_path = f"{data_dir}/1/masked_alignment_batch{batch_id}_{param_term}.maple"
    reset_file(maple_alignment_masked_file_path)
    path_store_n_masked = f"{data_dir}/1/n_masked_batch_{batch_id}_{param_term}.tsv"
    reset_file(path_store_n_masked)

    # Iterate over the test_list
    for path, sample_name in tqdm(sample_list,
                                desc="Processing samples",
                                total=n_samples,
                                unit="sample",
                                mininterval=90):
        clean_depth_list = []
        cons_depth_list = []
        ref_pos_list = []
        unma_seq = []
        n_het_sites = 0
        unmasked_pos_dict = build_masked_pos_dict() # This dict tracks all the positions masked by Viridian and our initial masking
        
        # unzip and open the file
        with smart_open(path, "rt") as f:
            # Read the lines and create the needed values to apply the thresholds
            # List of all the clean_depth, list of het_prop, median_depth
            
            # b'Ref_pos\tRef_nt\tCons_pos\tCons_nt\tMasked_cons_nt\tAmplicon\tPrimer\tMask\tTotal_depth\tClean_depth\tCons_depth\tA\ta\tC\tc\tG\tg\tT\tt\tI\ti\tD\td\tX_A\tX_a\tX_C\tX_c\tX_G\tX_g\tX_T\tX_t\tX_I\tX_i\tX_D\tX_d\n'
            
            # Iterate over the lines in the file
            next(f)  # Skip header line
            for line in f:
                # Decode the line and split it into columns
                line = line.strip().split("\t")
                # Extract the values from the columns
                if line[9].isdigit(): # clean_depth
                    # Extract the values from the columns
                    ref_pos = int(line[0]) 
                    clean_depth = int(line[9])
                    cons_depth = int(line[10])
                    # Lowercase the reference and unmasked sequences
                    cons_nt = line[4].lower() # Consensus is Viridian's masked consensus (Masked_cons_nt)
                    
                    # minor allele calculation
                    counts = [
                        int(line[11]) + int(line[12]),  # A
                        int(line[13]) + int(line[14]),  # C
                        int(line[15]) + int(line[16]),  # G
                        int(line[17]) + int(line[18])   # T
                    ]
                    if clean_depth > 50:                    
                        # find minor allele (second highest)
                        max_idx = counts.index(max(counts))
                        counts[max_idx] = 0
                        minor_idx = counts.index(max(counts))
                        minor_allele = "ACGT"[minor_idx]
                        het_depth = counts[minor_idx] / clean_depth
                        
                        if het_depth >= het_thr:
                            # Considered as a heterozygous site if there's a minor allele with proportion higher than het_thr
                            n_het_sites += 1

                else:
                    ref_pos = int(line[0])
                    clean_depth = 0
                    cons_depth = 0
                    cons_nt = 'n'
                    # Skip the first and last lines where depths col are filled with .
                
                if (len(ref_pos_list) > 0) and (ref_pos == ref_pos_list[-1]):
                    # Insertion
                    unma_seq[-1] += cons_nt
                
                else:
                    unma_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    clean_depth_list.append(clean_depth)
                    cons_depth_list.append(cons_depth)
                    
                    if cons_nt == 'n' and not unmasked_pos_dict[ref_pos]:
                        unmasked_pos_dict[ref_pos] = True # Ref_pos is 1-based so no need to do +1
        
        if len(unma_seq) != 29903:
            # MAPLE needs all the sequences to be of the same length
            # With the insertion management, normally this if condition should never be reached
            print(f"Warning: Sequence length is not 29903 for {sample_name} with len {len(unma_seq)}. Skipping this sample.")
            continue
        
        # MASKING
        # Compute median depth
        median_depth = median(clean_depth_list)
        # Compute dropout masking threshold for the sample
        depth_int_thr = median_depth * depth_thr
        # Retrieve the indexes where the het_prop is under the threshold AND the clean depth is under the threshold
        # Build the masked sequence and random masked sequence        
        unma_seq = list(unma_seq)
        # Masking of the very heterozygous or prone to errors positions according to Nicola's papers 
        unma_seq = apply_maple_masking(unma_seq, ref) # WE INPUT THE REFERENCE TO REPLACE WITH THE REF NT INSTEAD OF 'n'
        
        # PHYCD masking
        masked_seq = unma_seq.copy()
        masked_pos_dict = unmasked_pos_dict.copy() # This dict tracks all masked positions in the dropout masked sequence

        random_seq = unma_seq.copy()
        random_pos_dict = unmasked_pos_dict.copy() # This dict tracks all masked positions in the randomly masked sequence
        shift = random.randint(1000, 29000)
        
        for i, (clean_depth, cons_depth) in enumerate(zip(clean_depth_list, cons_depth_list)):
            if clean_depth <= max(depth_int_thr, 50):
                # 50 is the minimum depth for consensus calling in our pipeline
                if not masked_pos_dict[i+1]: # the dict are 1-based
                    # We only count the positions that are masked by DPCA and NOT MAPLE or VIRIDIAN
                    masked_seq[i] = "n"
                    masked_pos_dict[i+1] = True
                    if not random_pos_dict[((i + shift) % 29903) + 1]:
                        random_seq[((i + shift) % 29903)] = "n"
                        random_pos_dict[((i + shift) % 29903) + 1] = True
        
        # Check the number of masked positions and apply the threshold
        n_masked_masked = sum(1 for pos_masked in masked_pos_dict.values() if pos_masked)
        n_random_masked = sum(1 for pos_masked in random_pos_dict.values() if pos_masked)
        n_masked_cons = sum(1 for pos_masked in unmasked_pos_dict.values() if pos_masked)
        n_het_sites = int(n_het_sites)
        n_dropout_masked = n_masked_masked - n_masked_cons
        # Store these values in a tsv file
        with open(path_store_n_masked, "a") as f:
            f.write(f"{sample_name}\t{n_masked_cons}\t{n_masked_masked}\t{n_het_sites}\n")

        # KEEPING SAMPLES TO BUILD THE CLEAN TREE
        if n_het_sites <= max_n_het_sites and n_dropout_masked <= max_dropout_masked:
            # print(f"Warning: No masked positions for {sample_name}. Storing this sample for the clean tree.")
            
            with open(clean_tree_files_path, "a") as f:
                # SHOULD I USE unma_seq OR ORIGINAL SEQ? unma_seq has the MAPLE masking while original_seq is the exact Viridian consensus
                cons_entry = build_maple_entry(unma_seq, ref, f"{sample_name}")
                f.write(f"{cons_entry}\n")
                
            # This sample cannot be tested for contamination as it is used to build the clean tree
            continue
        
        # MAKING SURE THAT DROPOUT MASKED AND RANDOMLY MASKED SEQUENCES HAVE THE SAME NUMBER OF MASKED POSITIONS
        n_diff = n_masked_masked - n_random_masked
        # We want n_diff to be 0 after this loop
        while n_diff != 0:
            # Mask n_diff consecutive positions in the random sequence, if we fall on already masked positions, do it again with the new n_diff value
            start = random.randint(73, 29700) # 29769 and above are masked, 1 - 73 are masked
            for i in range(n_diff):
                if not random_pos_dict[((start + i) % 29903) + 1]:
                    random_seq[(start + i) % 29903] = "n"
                    random_pos_dict[((start + i) % 29903) + 1] = True
                    n_random_masked += 1
            
            n_diff = n_masked_masked - n_random_masked
        
        # Add the VDN unmasked consensus, our masked sequence and the random masked sequence to the output file
        with open(path_write_file, "a") as f:
            cons_entry = build_maple_entry(unma_seq, ref, f"unmasked_{sample_name}")
            masked_entry = build_maple_entry(masked_seq, ref, f"masked_{sample_name}")
            random_entry = build_maple_entry(random_seq, ref, f"random_{sample_name}")
            f.write(f"{cons_entry}\n{masked_entry}\n{random_entry}\n")
            
        with open(maple_alignment_unmasked_file_path, "a") as f:
            f.write(f"{cons_entry}\n")
        with open(maple_alignment_random_file_path, "a") as f:
            f.write(f"{random_entry}\n")
        with open(maple_alignment_masked_file_path, "a") as f:
            f.write(f"{masked_entry}\n")
        
        
        # Store the mutations masked in the masked and random sequences
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
     
        # Write the mutations masked to a file
        with open(path_mut_file, "a") as f:
            masked_mut = ';'.join(masked_mut)
            masked_mut = masked_mut.upper()
            
            rand_masked_mut = ';'.join(rand_masked_mut)
            rand_masked_mut = rand_masked_mut.upper()
            f.write(f"{sample_name}\t{masked_mut}\t{rand_masked_mut}\n")

    # Compress the final file if wanted
    if compress:
        compress_file(path_write_file, threads=1)

if __name__ == "__main__":
    # Generate the batches
    print(platform.python_implementation(), platform.python_version())
    print(f"Data directory: {data_dir}", flush=True)    
    if batch_id==0:
        # Remove the files in the batches folder if they exist
        batches_folder = Path(f"{data_dir}/1/")
        batches_folder.mkdir(parents=True, exist_ok=True)
        for file in batches_folder.glob(f"*"):
            file.unlink()
        print("Cleared the batches folder.")
    else:
        print("Waiting for process 0 to clear the batches folder...", flush=True)
        time.sleep(90)

    batchs_samples_list = generate_batchs_pairs(n_batch)
    sample_list = batchs_samples_list[batch_id]
    del batchs_samples_list
    
    print(param_term)
    # Run the function on each batch
    apply_masking(sample_list=sample_list,
                  het_thr=het_thr,
                  depth_thr=depth_thr,
                  max_dropout_masked=max_dropout_masked,
                  max_n_het_sites=max_n_het_sites,
                  path_ref_seq=path_ref_seq,
                  param_term=param_term,
                  data_dir=data_dir,
                  compress=compress,
                  )
    # Save the result
    print(f"Batch {batch_id} done")
