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
import pickle

parser = argparse.ArgumentParser(description='Detect potential contaminated areas depending on several thresholds applied on the Viridian samples.')

parser.add_argument('--n_batch', type=int, default=256,
                    help='Number of batches/jobs. Default is 256.')
parser.add_argument('--batch_id', type=int, default=0,
                    help='Batch id to process. Default is 0.')

parser.add_argument('--n_het_thr', type=int, default=3,
                    help='Threshold for the number of heterozygous positions with similar het prop for the samples to be labeled as contaminated. Default is 3.')

parser.add_argument('--path_ref_seq', type=str, default="/nfs/research/goldman/anoufa/data/NC_045512.2.fasta",
                    help='Path to the reference sequence. Default is /nfs/research/goldman/anoufa/data/NC_045512.2.fasta.')
args = parser.parse_args()

# Get the arguments
n_batch = args.n_batch
batch_id = args.batch_id
n_het_thr = args.n_het_thr
path_ref_seq = args.path_ref_seq

def generate_batchs_pairs(n_batchs):
    """Generate a 2D list with structure such that L[i][j] = sample j of batch i.
    Samples are pairs (path, read_name) where path is the path to the qc file.

    Args:
        n_batchs (int): Number of batches to generate.
    """
    path_vdn = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"
    # Generate a list of ids to chose samples
    batchs_samples_list = [[] for _ in range(n_batchs)]

    with lzma.open(path_vdn, "rt") as f:
        # Each line is composed of read_name and path to qc file
        # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
        for i, line in enumerate(f):
            if i == 0:
                continue
            # Get the read_name and path
            read_name, path = line.strip().split("\t")
            path = path + "/qc.tsv.gz"
            batchs_samples_list[i % n_batchs].append((path, read_name))
    
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
    maskedPoss[27507]=True
    maskedPoss[8835]=True # high het
    maskedPoss[15521]=True # high het
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
    maskedPoss[23118]=True # high het
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

def bin_het_prop(het_prop, het_prop_list):
    # Bin the heterozygosity proportion into 10 categories.
    # Categories should overlap to avoid edge cases.
    # We bin proportion of the second most abundant nucleotide
    
    if 0.06 <= het_prop < 0.10:
        het_prop_list[0] += 1
    if 0.10 <= het_prop < 0.15:
        het_prop_list[1] += 1
    if 0.15 <= het_prop < 0.20:
        het_prop_list[2] += 1
    if 0.20 <= het_prop < 0.25:
        het_prop_list[3] += 1
    if 0.25 <= het_prop < 0.30:
        het_prop_list[4] += 1
    if 0.30 <= het_prop < 0.35:
        het_prop_list[5] += 1
    if 0.35 <= het_prop < 0.40:
        het_prop_list[6] += 1
    if 0.40 <= het_prop < 0.45:
        het_prop_list[7] += 1
    if 0.45 <= het_prop < 0.50:
        het_prop_list[8] += 1
    return het_prop_list

def sort_samples(sample_list,
                path_write_file):
    """Read the Viridian qc file, apply several masking thresholds and
    generate a MAPLE file with the VDN consensus sequences, our masked sequences, and random masked sequences.

    Args:
        sample_list (list): list of tuples containing the sample name and the path to the file.
        het_thr (float): threshold for the heterozygous proportion.
        depth_thr (float): threshold for the depth.
        prop_under_depth_thr (float): threshold for the maximal proportion of positions with clean depth below the threshold.
        typical_depth_thr (int): median depth threshold.
        path_ref_seq (str): path to the reference sequence.
        path_write_file (str): path to the output file.
    Returns:
        str: path to the MAPLE alignment file.
    """
    
    n_samples = len(sample_list)
    masked_pos_dict = build_masked_pos_dict()
    # Iterate over the test_list
    for path, read_name in tqdm(sample_list,
                                desc="Processing samples",
                                total=n_samples,
                                unit="sample",
                                mininterval=300):
        clean_depth_list = []
        cons_depth_list = []
        ref_pos_list = []
        het_prop_list = [0] * 9
        cons_seq = []
        
        # unzip and open the file
        with gzip.open(path, "rb") as f:
            # Read the lines and create the needed values to apply the thresholds
            # List of all the clean_depth, list of het_prop, median_depth
            
            # b'Ref_pos\tRef_nt\tCons_pos\tCons_nt\tMasked_cons_nt\tAmplicon\tPrimer\tMask\tTotal_depth\tClean_depth\tCons_depth\tA\ta\tC\tc\tG\tg\tT\tt\tI\ti\tD\td\tX_A\tX_a\tX_C\tX_c\tX_G\tX_g\tX_T\tX_t\tX_I\tX_i\tX_D\tX_d\n'
            
            # Iterate over the lines in the file
            for line in f:
                # ignore header line
                if line.startswith(b"Ref_pos"):
                    continue
                # Decode the line and split it into columns
                line = line.decode("utf-8").strip().split("\t")
                # Extract the values from the columns
                try:
                    # Extract the values from the columns
                    ref_pos = int(line[0])
                    clean_depth = int(line[9])
                    cons_depth = int(line[10])
                    # Lowercase the reference and consensus sequences
                    cons_nt = line[4].lower()
                    
                    A_count = int(line[11]) + int(line[12])
                    C_count = int(line[13]) + int(line[14])
                    G_count = int(line[15]) + int(line[16])
                    T_count = int(line[17]) + int(line[18])
                    I_count = int(line[19]) + int(line[20])
                    D_count = int(line[21]) + int(line[22])
                    
                    
                except:
                    ref_pos = int(line[0])
                    clean_depth = 0
                    cons_depth = 0
                    cons_nt = line[4].lower()
                    
                    A_count, C_count, G_count, T_count = 0, 0, 0, 0
                    I_count, D_count = 0, 0
                    # Skip the first and last lines where depths col are filled with .
                
                if (len(ref_pos_list) > 0) and (ref_pos == ref_pos_list[-1]):
                    # Insertion
                    cons_seq[-1] += cons_nt
                
                else:
                    cons_seq.append(cons_nt)
                    ref_pos_list.append(ref_pos)
                    clean_depth_list.append(clean_depth)
                    cons_depth_list.append(cons_depth)
                    
                # Heterozygosity proportion (prop of the second most abundant nucleotide)
                total_count = A_count + C_count + G_count + T_count + I_count + D_count
                if total_count == 0:
                    # If no counts, skip this position
                    continue
                # Get the second most abundant nucleotide
                counts = [A_count, C_count, G_count, T_count, I_count, D_count]
                counts.sort(reverse=True)
                second_most_abundant = counts[1]
                if second_most_abundant == 0:
                    # If 100% consensus, skip this position
                    continue                
                
                het_prop = second_most_abundant / clean_depth
                
                if not masked_pos_dict[ref_pos] and cons_nt != 'N':
                    het_prop_list = bin_het_prop(het_prop, het_prop_list)
                    
        
        median_depth = median(clean_depth_list)
        
        if median_depth < 100:
            # If median depth is below 100, skip this sample
            continue
        
        
        # Find if the sample has a certain het_prop category coming back very often
        het_prop_labels = ['0.06-0.10', '0.10-0.16', '0.15-0.21', '0.20-0.26',
                          '0.25-0.31', '0.30-0.36', '0.35-0.41', '0.40-0.46', '0.45-0.51']
        
        het_prop_max = max(het_prop_list)
        het_prop_idx = het_prop_list.index(het_prop_max)
        het_prop_label = het_prop_labels[het_prop_idx]
        
        if het_prop_max > n_het_thr:
            # Check if the other proportions are below the threshold
            if all(x <= n_het_thr for i, x in enumerate(het_prop_list) if i != het_prop_idx):
                # If the maximum proportion is above the threshold and the others are below, consider this sample as potentially contaminated
                with open(path_write_file, "a") as f:
                    f.write(f"{read_name}\t{het_prop_label}\t{het_prop_max}\t{het_prop_list}\n")
    
    
    return path_write_file

if __name__ == "__main__":
    # Generate the batchs
    print(platform.python_implementation(), platform.python_version())

    batchs_samples_list = generate_batchs_pairs(n_batch)
    
    sample_list = batchs_samples_list[batch_id]
    
    del batchs_samples_list
    
    path_write_file = f"/nfs/research/goldman/anoufa/data/dpca/batches/sort_samples_{batch_id}_{n_het_thr}.txt"
    
    # Run the function on each batch
    path_write_file = sort_samples(sample_list,
                                    path_write_file)
    
    # Save the result
    print(f"Batch {batch_id} done")
    
