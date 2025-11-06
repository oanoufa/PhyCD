from pathlib import Path
import shutil
# import tqdm

# Linux paths (adjust as needed)
print("Current working directory:", Path.cwd())
source_dir = Path("/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Reads/")  # Base directory containing RR** folders
dest_dir = Path("/homes/anoufa/VDN_READS")  # Where to save renamed files
dest_dir.mkdir(parents=True, exist_ok=True)  # Create destination with parents

# Process all matching .gz files
for top_folder in ['D']:
    files = source_dir.glob(f"{top_folder}/RR[0-9]*/[0-9][0-9]/[0-9][0-9]/vdn.v1.0.0/*.gz")
    for gz_file in files:
        
        if gz_file.suffixes == ['.json', '.gz']:
            continue
        
        # Extract path components
        parts = gz_file.parts
        
        # Construct ID (e.g., "DRR010203" from "D/RR01/02/03/...")
        top_folder = gz_file.parts[-6]  # Gets "D", "E", or "S"
        rr_num = parts[-5][2:]  # Gets "****" from "RR****"
        sub1 = parts[-4]        # Gets first ** ("02")
        sub2 = parts[-3]        # Gets second ** ("03")
        file_id = f"{top_folder}RR{rr_num}{sub1}{sub2}"
        
        # Create new filename and copy
        new_name = f"{file_id}_{gz_file.name}"
        shutil.copy2(gz_file, dest_dir / new_name)
        
        # print(f"Processed: {gz_file} → {new_name}")

print(f"Finished! Files saved to {dest_dir}")