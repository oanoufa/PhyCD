# Auxiliary functions used in several scripts
from pathlib import Path
import shutil
import gzip
import lzma
import _params
from tqdm import tqdm


def generate_sh_param_file():
    """Generate a sh param file containing the parameters in the python params.py file.
    """
    with open("/nfs/research/goldman/anoufa/shell_scripts/params.sh", "w") as f:
        for key, value in _params.__dict__.items():
            if not key.startswith("_") and not callable(value):
                print(f"Writing {key}={value} to params.sh")
                f.write(f"{key}={value}\n")

def update_params_file(params_file: str, updates: dict):
    """
    Update or add parameters in a Python params file.

    Args:
        params_file (str): Path to the .py params file (e.g., "_params.py").
        updates (dict): A dictionary of {parameter_name: value} pairs.
                        Values will be written as Python literals.
    """
    path = Path(params_file)

    if path.exists():
        with open(path, "r") as f:
            lines = f.readlines()
    else:
        lines = []

    updated_keys = set()

    for i, line in enumerate(lines):
        stripped = line.strip()
        for key, val in updates.items():
            if stripped.startswith(f"{key} "):  # e.g., "processed_placement = ..."
                lines[i] = f"{key} = {repr(val)}\n"
                updated_keys.add(key)
                break

    for key, val in updates.items():
        if key not in updated_keys:
            lines.append(f"{key} = {repr(val)}\n")

    with open(path, "w") as f:
        f.writelines(lines)
        
def compress_file(path):
    path = Path(path)
    # Compress the file with lzma or gzip
    print(f"Compressing {path}...")
    with open(path, "rb") as f_in:
        gz_path = path.parent / (path.name + ".gz")
        with gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    # Remove the original file
    path.unlink()
    
# Build MAPLE format entry for the unmasked and masked sequences
def build_maple_entry(seq, ref, seq_name):
    """Code adapted from Nicola De Maio
    This code builds a MAPLE format entry to insert in MAPLE file for a given sequence and reference.
    """
    
    # Assertion errors to check we are fed lists
    assert isinstance(seq, list), "seq should be a list"
    entry = [f">{seq_name}"]
    state = 0
    length = 0
    lRef = len(ref)
    seq_list = []
    
    for i in range(lRef):
        current_nuc = seq[i]
        
        if len(current_nuc) > 1:
            # Insertion --> replace with the first nucleotide
            # print(f"Warning: Insertion {current_nuc} masked at position {i+1} in {seq_name}.")
            current_nuc = current_nuc[0]

        if state == 1:
            # State 1: n
            if current_nuc == "n":
                length += 1
            else:
                # Close the n run
                seq_list.append(("n", i + 1 - length, length))
                length = 0
                state = 0
        elif state == 2:
            # State 2: -
            if current_nuc == "-":
                length += 1
            else:
                # Close the - run
                seq_list.append(("-", i + 1 - length, length))
                length = 0
                state = 0
        if current_nuc == "n" and state != 1:
            # Start a new n run
            length = 1
            state = 1
        elif current_nuc == "-" and state != 2:
            # Start a new - run
            length = 1
            state = 2
        elif current_nuc != ref[i] and current_nuc != "-" and current_nuc != "n":
            # Write the differing nucleotide
            seq_list.append((current_nuc, i + 1))

    # Close any run at the end
    if state == 1:
        seq_list.append(("n", lRef + 1 - length, length))
    elif state == 2:
        seq_list.append(("-", lRef + 1 - length, length))

    for m in seq_list:
        if len(m) == 2:
            entry.append(f"{m[0]}\t{m[1]}")
        else:
            entry.append(f"{m[0]}\t{m[1]}\t{m[2]}")
    return "\n".join(entry)

def change_ref_maple(path_mpl, new_ref_seq):
    """Change the reference sequence in a MAPLE file."""
    
    new_mpl_path = path_mpl.replace(".maple", "_NC_045512.2.maple")
    with open(new_mpl_path, "w") as f2:
        f2.write(f">REF\n{new_ref_seq}\n")
        
    with open(path_mpl, "r") as f:
        f.readline()  # Skip the first line (reference sequence)
        former_ref = f.readline().strip()
        start = False
        for line in tqdm(f, desc="Processing MAPLE file reference modification", mininterval=30):
            if line.startswith(">"):
                if start:
                    # Save the previous sequence
                    with open(new_mpl_path, "a") as f2:
                        maple_entry = build_maple_entry(current_seq, new_ref_seq, seq_name[1:])
                        f2.write(maple_entry + "\n")
                        
                # Start a new sequence
                start = True
                current_seq = list(former_ref)  # Start with the former reference sequence
                seq_name = line.strip()
            
            else:
                if line.startswith("N") or line.startswith("-"):
                    parts = line.split("\t")
                    nuc = parts[0]
                    pos = int(parts[1]) - 1
                    num = parts[2]
                    
                    for i in range(int(num)):
                        current_seq[pos + i] = nuc
                    
                if line.startswith("A") or line.startswith("C") or line.startswith("G") or line.startswith("T"):
                    parts = line.split("\t")
                    nuc = parts[0]
                    pos = int(parts[1]) - 1
                    current_seq[pos] = nuc
                    
            
        # Save the last sequence
        if start:
            with open(new_mpl_path, "a") as f2:
                maple_entry = build_maple_entry(current_seq, new_ref_seq, seq_name[1:])
                f2.write(maple_entry + "\n")
                
    return new_mpl_path


def build_maple_file(path_ref_seq, path_write_file):
    
    shutil.copy(path_ref_seq, path_write_file)
    
    with open(path_ref_seq, "r") as f:
        # Reference sequence is at line 2
        f.readline()
        ref_seq = f.readline().strip()
    
    return ref_seq

def parse_maple_file(path_mpl):
    """Iterate through a MAPLE file and output the reference and a dict that links the entry name to its entry.
    The dict has the following structure:
    {
        "seq_name1": [entry_lines],
        "seq_name2": [entry_lines],
        ...
    }

    Args:
        path_mpl (str): Path to the MAPLE file.
        
    Returns:
        ref_name (str): Reference sequence name.
        ref_seq (str): Reference sequence.
        maple_dict (dict): Dictionary mapping sequence names to their entries.
    """
    
    maple_dict = {}
    with open(path_mpl, "r") as f:
        ref_name = f.readline().strip()[1:]  # Skip the '>' character
        ref_seq = f.readline().strip()
        
        current_seq_name = None
        current_entry = []
        
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save the previous entry
                if current_seq_name is not None:
                    maple_dict[current_seq_name] = current_entry
                
                # Start a new entry
                current_seq_name = line[1:]  # Skip the '>' character
                current_entry = []
            else:
                current_entry.append(line)
        
        # Save the last entry
        if current_seq_name is not None:
            maple_dict[current_seq_name] = current_entry
            
    return ref_name, ref_seq, maple_dict

def expand_ambiguous_mutation(mut, remove_starting_nt = False):
    """
    Expand an ambiguous mutation like C18745W to all possible unambiguous mutations (here C18745A and C18745T).
    AMBIGUOUS SIGNS: W (A or T), S (C or G), M (A or C), K (G or T), R (A or G), Y (C or T)
    
    input: mut (str): mutation string 'C18745W'
    output: list of expanded mutations (list of str): ['C18745A', 'C18745T']
    """
    ambiguous_dict = {
        "W": ["A", "T"],
        "S": ["C", "G"],
        "M": ["A", "C"],
        "K": ["G", "T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
    }

    if len(mut) < 2:
        return [mut]
    
    if remove_starting_nt:
        mut = mut[1:]

    base = mut[:-1]
    ambiguous_base = mut[-1]

    if ambiguous_base in ambiguous_dict:
        return [base + unambiguous for unambiguous in ambiguous_dict[ambiguous_base]]

    return [mut]

def generate_sample_list(sample_names):
    """Get the sample names in the result dataset and retrieve the paths to their qc.tsv.gz files.
    """
    path_vdn = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"
    # if sample_names is a list, convert to set for faster lookup
    if isinstance(sample_names, list):
        sample_names = set(sample_names)
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

