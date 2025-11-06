# Generate distributions over the whole dataset

import platform
from tqdm import tqdm
import gzip
from statistics import median
import lzma
import argparse
import shutil
from pathlib import Path
import pickle

parser = argparse.ArgumentParser(description='Detect potential contaminated areas depending on several thresholds applied on the Viridian samples.')

parser.add_argument('--n_batch', type=int, default=256,
                    help='Number of batches/jobs. Default is 256.')
parser.add_argument('--batch_id', type=int, default=0,
                    help='Batch id to process. Default is 0.')

parser.add_argument('--het_thrs', type=float, nargs='+', default=[0.95, 0.9, 0.8, 0.7],
                    help='Thresholds for the heterozygous proportion. Default is [0.95, 0.9, 0.8, 0.7].')
parser.add_argument('--depth_thrs', type=float, nargs='+', default=[0.2, 0.1, 0.05, 0.02, 0.01],
                    help='Thresholds for the clean depth. Default is [0.1, 0.02].')
args = parser.parse_args()

# Get the arguments
n_batch = args.n_batch
batch_id = args.batch_id
het_thrs = args.het_thrs
depth_thrs = args.depth_thrs

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

def bin_length_values(length):
    # Bin consensus length values with cat for 29903 to 29915
    cat = length-29903
    if cat > 12:
        cat = 12
    return int(cat)

def bin_median_values(median):
    # Bin values with cat for each hundred
    cat = median // 100
    
    if cat > 12:
        cat = 12
        
    return int(cat)

def bin_n_het_sites(n_het_sites):
    # Bin values with one cat until 5 then 5 - 9, 10 - 19, 20 - 29, 30+
    if n_het_sites < 5:
        cat = n_het_sites
    elif n_het_sites < 10:
        cat = 5
    elif n_het_sites < 20:
        cat = 6
    elif n_het_sites < 30:
        cat = 7
    else:
        cat = 8
    
    return int(cat)

def bin_prop_below_depth(prop_below_depth):
    # Bin values with cat for each 0.02
    cat = prop_below_depth*100 // 2
    
    if cat >= 20:
        cat = 19
        
    return int(cat)

def bin_amplicon_length(length):
    # Bin values with cat for each 100 starting at 50
    
    cat = (length - 50) // 100
    return int(cat)

def compute_dropout_lengths(dropout_list):
    lengths = []
    count = 0
    for state in dropout_list:
        if state:
            count += 1
        elif count > 0:
            lengths.append(count)
            count = 0
    if count > 0:  # Handle case where list ends with True
        lengths.append(count)
    return lengths

def bin_dropout_length(length):
    # Bin into categories every 50
    cat = length // 50
    return int(cat)

def gen_metrics(sample_list,
                het_thrs,
                depth_thrs):
    """Generate counts for different metrics.

    Args:
        sample_list (list): list of tuples containing the sample name and the path to the file.
        het_thrs (list): list of thresholds for the heterozygous proportion.
        depth_thrs (float, optional): _description_. Defaults to 0.2.
    """
    
    n_samples = len(sample_list)
    median_depth_list = [0] * 13
    mean_depth_list = [0] * 13
    n_het_thr = len(het_thrs)
    het_sites_list = [[0] * 9 for _ in range(n_het_thr)]
    n_depth_thr = len(depth_thrs)
    prop_below_depth_list = [[0] * 20 for _ in range(n_depth_thr)]
    length_list = [0] * 13
    amplicons_length_list = [0] * 20
    dropout_lengths_list = [0] * 400
    # Iterate over the test_list
    for path, _ in tqdm(sample_list,
                                desc="Processing samples",
                                total=n_samples,
                                unit="sample",
                                mininterval=300):
        clean_depth_list = []
        cons_depth_list = []
        cons_seq = ""
        amplicons_dict = {}
        dropout = False
        dropout_states_list = []
        
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
                    clean_depth = int(line[9])
                    cons_depth = int(line[10])
                    # Lowercase the reference and consensus sequences
                    cons_seq += line[4].lower()
                except:
                    clean_depth = 0
                    cons_depth = 0
                    cons_seq += line[4].lower()
                    # Skip the first and last lines where depths col are filled with .
                clean_depth_list.append(clean_depth)
                cons_depth_list.append(cons_depth)
                
                amplicons = line[5]
                
                for amplicon in amplicons.split(";"):
                    if amplicon not in amplicons_dict:
                        amplicons_dict[amplicon] = 1
                    else:
                        amplicons_dict[amplicon] += 1
                
                if 0 < clean_depth <= 30:
                    # Enter dropout state if clean depth is below 30 (VDN threshold not to call consensus)
                    # AND if clean depth is not 0 (we want to check dropouts with a bit of depth that could be a contamination)
                    # Only leave dropout state if clean depth is above 50
                    if not dropout:
                        dropout = True
                elif clean_depth > 50:
                    dropout = False
                    
                dropout_states_list.append(dropout)
        
        cons_seq_len = len(cons_seq)
        length_list[bin_length_values(cons_seq_len)] += 1
        
        median_depth = median(clean_depth_list)
        mean_depth = sum(clean_depth_list) / len(clean_depth_list)
        # Store the mean depth
        mean_depth_list[bin_median_values(mean_depth)] += 1
        # Store the median depth
        median_depth_list[bin_median_values(median_depth)] += 1
        
        # Store the amplicon length
        amplicon_average_length = sum(amplicons_dict.values()) / len(amplicons_dict)
        amplicons_length_list[bin_amplicon_length(amplicon_average_length)] += 1
        
        # Store the dropout lengths
        dropout_lengths = compute_dropout_lengths(dropout_states_list)
        # Remove first and last dropout lengths (beginning and end of the sequence)
        if len(dropout_lengths) > 2:
            dropout_lengths = dropout_lengths[1:-1]
            for dropout_length in dropout_lengths:
                dropout_lengths_list[bin_dropout_length(dropout_length)] += 1
        
    
        # 2. Check if the proportion of clean depth below the threshold is above the threshold
        for i, depth_thr in enumerate(depth_thrs):
            depth_int_thr = median_depth * depth_thr
            prop_under_depth = sum(1 for d in clean_depth_list if d < depth_int_thr) / len(clean_depth_list)
        
            prop_below_depth_list[i][bin_prop_below_depth(prop_under_depth)] += 1
        # 3. Retrieve the indexes where the het_prop is under the threshold AND the clean depth is under the threshold
        # Build the masked sequence
        n_het_sites = [0] * n_het_thr

        for _, (clean_depth, cons_depth) in enumerate(zip(clean_depth_list, cons_depth_list)):
            het_prop = 1 if clean_depth == 0 else cons_depth / clean_depth
            for i in range(n_het_thr):
                if het_prop <= het_thrs[i]:
                    n_het_sites[i] += 1

        cat_n_het_sites = [bin_n_het_sites(n) for n in n_het_sites]
        for i in range(len(cat_n_het_sites)):
            # Store the number of het sites
            het_sites_list[i][cat_n_het_sites[i]] += 1
    
    storing_list = [median_depth_list, mean_depth_list, het_sites_list, prop_below_depth_list, length_list, amplicons_length_list, dropout_lengths_list]
                
    return storing_list


if __name__ == "__main__":
    # Generate the batchs
    print(platform.python_implementation(), platform.python_version())

    batchs_samples_list = generate_batchs_pairs(n_batch)
    
    sample_list = batchs_samples_list[batch_id]
    
    del batchs_samples_list
    
    print(het_thrs, depth_thrs)
        
    # Run the function on each batch
    storing_list = gen_metrics(sample_list,
                                het_thrs=het_thrs,
                                depth_thrs=depth_thrs)
    
                                                
    # Access the total list and add the results to the file
    path_storing_list = Path(f"/nfs/research/goldman/anoufa/data/dpca/gen_metrics/storing_file_{batch_id}.pkl")
    
    with open(path_storing_list, "wb") as f:
        pickle.dump(storing_list, f)
        print(f"Storing list saved to {path_storing_list}")
    
    # Save the result
    print(f"Batch {batch_id} done")
    
