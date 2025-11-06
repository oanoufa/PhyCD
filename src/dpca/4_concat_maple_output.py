from pathlib import Path
import pandas as pd
from tqdm import tqdm
import _params
from _aux_functions import compress_file

n_batch = _params.n_batch
het_thr = _params.het_thr
depth_thr = _params.depth_thr
prop_under_depth_thr = _params.prop_under_depth_thr
typical_depth_thr = _params.typical_depth_thr
n_masked_thr = _params.n_masked_thr
path_ref_seq = _params.path_ref_seq
inputTree = _params.inputTree
inputRates = _params.inputRates
n_diff_mut = _params.n_diff_mut
masked_max_dist = _params.masked_max_dist
masking_ratio = _params.masking_ratio
param_term = _params.param_term


if __name__ == "__main__":
    

    # Access the total median cat list and add the results to the file
    final_path = Path(f"/nfs/research/goldman/anoufa/data/MAPLE_output/output_FULL_metaData_samplePlacements_{param_term}.tsv")

    final_mut_tsv = pd.DataFrame(columns=["sample", "placements", "optimizedBlengths", "mutations"])
    
    batches_folder_path = Path("/nfs/research/goldman/anoufa/data/MAPLE_output/batches/")
    
    dict_file_checker = {i: False for i in range(n_batch)}
    
    for mut_tsv_file in tqdm(batches_folder_path.glob("*_samplePlacements.tsv"),
                                desc="Processing files", unit="file",
                                mininterval=180):
        
        # Load the data
        mut_tsv = pd.read_csv(mut_tsv_file, sep="\t")
        # Add the data to the final dataframe
        final_mut_tsv = pd.concat([final_mut_tsv, mut_tsv], ignore_index=True)
        
        # Check if all n_batch files are processed
        # TYPICAL FILE NAME /nfs/research/goldman/anoufa/data/MAPLE_output/batches/output_500_metaData_samplePlacements.tsv
        batch_id = int(mut_tsv_file.stem.split("_")[1])
        dict_file_checker[batch_id] = True
        
        mut_tsv_file.unlink()
    
    if not all(dict_file_checker.values()):
        missing_batches = [k for k, v in dict_file_checker.items() if not v]
        print(f"Not all batch files were processed. Missing batches: {missing_batches}")

    else:
        print(f"All {n_batch} batch files processed successfully.")
    
    final_mut_tsv.to_csv(final_path, sep="\t", index=False)
    
    print("Concatenation complete. Now compressing the final file...", flush=True)
    
    # Compress the final file
    compress_file(final_path)
    
    print("Compression complete. Now deleting all updatedBlengths.tree files...", flush=True)
    
    # Delete all .tree files in the directory
    for tree_file in batches_folder_path.glob("*updatedBlengths.tree"):
        tree_file.unlink()
        
    print("All .tree files deleted.")
    # Delete the MAPLE alignment files if everything went well
    maple_files_path = Path("/nfs/research/goldman/anoufa/data/MAPLE_input/1_batches/")
    print("Now deleting all MAPLE alignment files...", flush=True)
    for maple_file in maple_files_path.glob("maple_alignment_batch*"):
        # Delete all the files except the ones in missing batches
        batch_id = int(maple_file.stem.split("batch")[1])
        if dict_file_checker.get(batch_id, True):
            maple_file.unlink()
    
    print(f"Final output saved to {final_path.with_suffix('.gz')}")
    