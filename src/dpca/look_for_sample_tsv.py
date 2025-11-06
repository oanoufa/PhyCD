from tqdm import tqdm
from pathlib import Path

if __name__ == "__main__":
    
    sample_to_look_for = "SRR19675317"
    sample_types_dict = {
        "consensus": False,
        "masked": False,
        "random": False
    }
    # Define the path to the folder containing the tsv files
    folder_path = Path("/nfs/research/goldman/anoufa/data/MAPLE_output/")

    mutations_dict = {}
    # Iterate over the tsv files in the folder
    for file_path in folder_path.glob("*_samplePlacements.tsv"):
        
        new_df_list = []
        # Read the tsv file
        with open(file_path, "r") as file:
        
            # Extract the sample lengths for the three types of sequences
            for line in file:
                # Skip first line (header) and empty lines
                line = line.strip()
                
                if not line or line.startswith("sample"):
                    continue
                # Columns:  index   sample	                placements	                    optimizedBlengths (<topBlength>/<bottomBlength>/<sampleBlength>)   mutations
                # Example line: 0	consensus_ERR4239172	in2495143:0.9997103398863393	in2495143:(6.68997107947098e-05/0/0)	T10029C;G21618C;G22917T;A22995C;T23063A;A23604...
                
                # We are interested in the sampleBlengths of the placement with highest support.
                columns = line.split("\t")
                sample = columns[0]
                try:
                    sample_type, sample_name = sample.split("_", 1)  # Split the sample name to get the type and read name
                except:
                    print(sample, columns, line, flush=True)
                    continue
                
                if sample_name != sample_to_look_for:
                    continue
                
                print(line)
                sample_types_dict[sample_type] = True
                
        if sample_types_dict["consensus"] and sample_types_dict["masked"] and sample_types_dict["random"]:
            print("Found all sample types for", sample_to_look_for)
            break