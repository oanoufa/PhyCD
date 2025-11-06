from tqdm import tqdm
from pathlib import Path
import gzip
import lzma
import _params

param_term = _params.param_term

if __name__ == "__main__":
    
    sample_to_look_for = "ERR9417311"
    sample_type = 'unmasked'

    print(f"Looking for sample {sample_to_look_for} in alignment file {sample_type}")
    # Iterate over the tsv files in the folder
    file_to_look_for = f"/nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/{sample_type}_alignment_{param_term}.maple"

    # Read the tsv file
    with open(file_to_look_for, "r") as file:
    # Read the tsv file
    
        # Extract the sample lengths for the three types of sequences
        for line in file:
            # Skip first line (header) and empty lines
            line = line.strip()

            if line == f">{sample_type}_{sample_to_look_for}":
                line = next(file).strip()
                while not line.startswith(">"):
                    print(line)
                    line = next(file).strip()
                print(f"Found {sample_type} entry for {sample_to_look_for} in {file_to_look_for}")
                break