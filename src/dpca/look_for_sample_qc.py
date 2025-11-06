# The goal of this script is to find a sample in the dataset and download his qc file back to our data folder

from pathlib import Path
import shutil
import lzma
import gzip


if __name__ == "__main__":
    path_vdn = Path("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz")
    
    sample_names = ["SRR17439727"]

    for sample_name in sample_names:
        path_to_copy_to = Path("/nfs/research/goldman/anoufa/data/samples/" + sample_name + "/")
        
        # Create the directory if it does not exist
        path_to_copy_to.mkdir(parents=True, exist_ok=True)

        with lzma.open(path_vdn, "rt") as f:
            # Each line is composed of read_name and path to qc file
            # Iterate over the lines and get the path, read_name if the line's index is in chosen_samples
            for i, line in enumerate(f):
                if i == 0:
                    continue
                # Get the read_name and path
                read_name, path = line.strip().split("\t")
                path = path + "/qc.tsv.gz"
                if read_name == sample_name:
                    # Copy the file to the data folder
                    path = Path(path)
                    path_end = sample_name + "_qc.tsv"
                    
                    path_to_copy_to = path_to_copy_to / path_end
                    if path.exists():
                        # Ungzip the file
                        print(f"Unzipping {path} to {path_to_copy_to}")
                        with gzip.open(path, "rb") as f_in:
                            with open(path_to_copy_to, "wb") as f_out:
                                shutil.copyfileobj(f_in, f_out)
                    else:
                        print(f"Path {path} does not exist")
                    break