# Auxiliary functions used in several scripts
from pathlib import Path
import shutil
import lzma
from tqdm import tqdm

def compress_file(path, level=3, threads=-1):
    """Compress a file using Zstandard (.zst)."""
    import zstandard as zstd
    from pathlib import Path
    path = Path(path)

    print(f"Compressing {path} with Zstandard...")

    cctx = zstd.ZstdCompressor(level=level, threads=threads)

    zst_path = path.parent / (path.name + ".zst")

    with open(path, "rb") as f_in, open(zst_path, "wb") as f_out:
        with cctx.stream_writer(f_out) as compressor:
            shutil.copyfileobj(f_in, compressor)
            compressor.flush(zstd.FLUSH_FRAME)

    # Optionally remove original file
    path.unlink()

    return zst_path

def smart_open(path, mode="rt"):
    """Open a file automatically handling gzip (.gz), zstd (.zst) and lzma (.xz) compressed files based on their extension.
    """
    from pathlib import Path

    path = Path(path)
    if not path.exists():
        matches = list(path.parent.glob(path.name + "*"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected exactly one file matching {path}*, found {len(matches)}"
            )
        path = matches[0]

    suffix = path.suffix
    # --- Zstandard (.zst) ---
    if suffix == ".zst":
        import zstandard as zstd
        return zstd.open(path, mode)
    # --- Gzip (.gz) ---
    if suffix == ".gz":
        import gzip
        return gzip.open(path, mode)
    # --- LZMA (.xz) ---
    if suffix == ".xz":
        import lzma
        return lzma.open(path, mode)
    # --- Plain file ---
    return open(path, mode)

    
# Build MAPLE format entry for the unmasked and masked sequences
def build_maple_entry(seq, ref, seq_name):
    """Code adapted from Nicola De Maio
    This code builds a MAPLE format entry to insert in a MAPLE file for a given sequence and reference.
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

def build_maple_file(path_ref_seq, path_write_file):
    """Takes as input the path to the reference sequence and the path to the file to generate.
    Copies the content of the reference sequence to the file, and output the reference sequence as a list.
    """
    
    shutil.copy(path_ref_seq, path_write_file)
    
    with open(path_ref_seq, "r") as f:
        # Reference sequence is at line 2
        f.readline()
        ref_seq = f.readline().strip()
    
    return ref_seq

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

def generate_sample_list(sample_names, samples_dir):
    """Get the sample names in the result dataset and retrieve the paths to their qc.tsv.gz files.
    """
    # if sample_names is a list, convert to set for faster lookup
    if isinstance(sample_names, list):
        sample_names = set(sample_names)
    # Generate a list of ids to chose samples
    batchs_samples_list = []

    with lzma.open(samples_dir, "rt") as f:
        # Each line is composed of sample_name and path to qc file
        # Iterate over the lines and get the path, sample_name if the line's index is in chosen_samples
        next(f)  # Skip header line
        for line in f:
            # Get the sample_name and path
            sample_name, path = line.strip().split("\t")
            path = path + "/qc.tsv.gz"
            if sample_name in sample_names:
                batchs_samples_list.append((path, sample_name))
                sample_names.remove(sample_name)
            if len(sample_names) == 0:
                break
    print(f"Generated sample list of {len(batchs_samples_list)} entries", flush=True)
    return batchs_samples_list

def save_pickle_dict(d, path, compress, level=3, threads=-1):
    """Save a dict as a pickle file and compress it if compress is True.
    """
    import pickle
    from pathlib import Path
    path = Path(path)
    if compress:
        import zstandard as zstd
        if path.suffix != ".zst":
            path = path.with_suffix(path.suffix + ".zst")

        cctx = zstd.ZstdCompressor(level=level, threads=threads)

        with open(path, "wb") as f:
            with cctx.stream_writer(f) as compressor:
                pickle.dump(d, compressor, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(path, "wb") as f:
            pickle.dump(d, f, protocol=pickle.HIGHEST_PROTOCOL)

    return path

def load_pickle_dict(path, compress):
    """Load pickle dict (potentially compressed)
    """
    from pathlib import Path
    import pickle
    path = Path(path)
    if compress:
        matches = list(path.parent.glob(path.name + "*"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected exactly one file matching {path.name}*, found {len(matches)}"
            )
        path = matches[0]
    with smart_open(path, mode="rb") as f:
        return pickle.load(f)