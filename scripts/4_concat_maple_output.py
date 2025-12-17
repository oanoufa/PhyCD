from pathlib import Path
from tqdm import tqdm
from _aux_functions import compress_file
import shutil
import argparse


parser = argparse.ArgumentParser(description='Process the output of the sample placement batched script.')
parser.add_argument('--n_batch', type=int,
                    help='Number of batches to split the samples into.')
parser.add_argument('--data_dir', type=str,
                    help='Directory containing all the files generated during the pipeline.')
parser.add_argument('--param_term', type=str,
                    help='str built from the initial parameters used to recognize the files.')
parser.add_argument('--compress', type=int,
                    help='1 or 0, enables compression.')
args = parser.parse_args()

# Get the arguments
n_batch = args.n_batch
data_dir = args.data_dir
param_term = args.param_term
compress = args.compress

if __name__ == "__main__":
    # CONCATENATE MAPLE SAMPLE PLACEMENTS
    batches_folder_path = Path(f"{data_dir}/3/")
    dict_file_checker = {i: False for i in range(n_batch)}

    final_path = Path(f"{data_dir}/4/output_FULL_metaData_samplePlacements_{param_term}.tsv")
    if not final_path.exists():
        with open(final_path, "w") as outfile:
            outfile.write("sample\tplacements\toptimizedBlengths\tmutations\n")
        with open(final_path, "a") as outfile:
            for placement_file in tqdm(batches_folder_path.glob("*_samplePlacements.tsv"),
                                        desc="Processing files", unit="file",
                                        mininterval=180):
                with open(placement_file, "r") as infile:
                    infile.readline()  # skip header
                    shutil.copyfileobj(infile, outfile)
                batch_id = int(placement_file.stem.split("_")[1])
                dict_file_checker[batch_id] = True
                placement_file.unlink()
    
    if not all(dict_file_checker.values()):
        missing_batches = [k for k, v in dict_file_checker.items() if not v]
        print(f"Not all batch files were processed. Missing batches: {missing_batches}")
    else:
        print(f"All {n_batch} batch files processed successfully.")
    
    print("Concatenation complete. Now compressing the final file...", flush=True)
    
    # Compress the final file
    if final_path.exists() and compress:
        compress_file(final_path)
    print("Compression complete. Now deleting all updatedBlengths.tree files...", flush=True)
    
    # Delete all .tree files in the directory
    for tree_file in batches_folder_path.glob("*updatedBlengths.tree"):
        tree_file.unlink()
    print("All .tree files deleted.")
    
    # Delete the MAPLE alignment files if everything went well
    maple_files_path = Path(f"{data_dir}/1/")
    print("Now deleting all MAPLE alignment files...", flush=True)
    for maple_file in maple_files_path.glob("maple_alignment_batch*"):
        # Delete all the files except the ones in missing batches
        batch_id = int((maple_file.stem.split("batch")[1]).split('_')[0])
        if dict_file_checker.get(batch_id, True):
            maple_file.unlink()
    
    print(f"Final output saved to {final_path}")
    