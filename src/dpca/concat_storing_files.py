import pickle
from pathlib import Path


if __name__ == "__main__":
    # Generate the batchs
                                                
    # Access the total median cat list and add the results to the file
    final_path = Path("/nfs/research/goldman/anoufa/data/dpca/storing_file.pkl")
    
    folder_path = Path("/nfs/research/goldman/anoufa/data/dpca/gen_metrics/")
    

    # Generate a new file to reset it
    median_depth_list_total = [0] * 13
    mean_depth_list_total = [0] * 13
    n_het_sites_thrs = 4
    het_sites_list_total = [[0] * 9 for _ in range(n_het_sites_thrs)]
    n_depth_thrs = 5
    prop_below_depth_list_total = [[0] * 20 for _ in range(n_depth_thrs)]
    length_list_total = [0] * 13
    amplicons_length_list_total = [0] * 20
    dropout_lengths_list_total = [0] * 400
    
    for pkl_file in folder_path.glob("*.pkl"):
        batch_id = pkl_file.stem.split("_")[-1]
        print(f"Processing batch {batch_id}")
        
        path = folder_path / f"storing_file_{batch_id}.pkl"
        
        # Load the data
        with open(path, "rb") as f:
            storing_list = pickle.load(f)
            
        median_depth_list = storing_list[0]
        mean_depth_list = storing_list[1]
        het_sites_list = storing_list[2]
        prop_below_depth_list = storing_list[3]
        length_list = storing_list[4]
        amplicons_length_list = storing_list[5]
        dropout_lengths_list = storing_list[6]
        
        # Remove the files
        path.unlink()
        
        # Add the results to the total lists
        median_depth_list_total = [x + y for x, y in zip(median_depth_list_total, median_depth_list)]
        mean_depth_list_total = [x + y for x, y in zip(mean_depth_list_total, mean_depth_list)]
        
        length_list_total = [x + y for x, y in zip(length_list_total, length_list)]
        for i in range(n_het_sites_thrs):
            het_sites_list_total[i] = [x + y for x, y in zip(het_sites_list_total[i], het_sites_list[i])]

        for i in range(n_depth_thrs):
            prop_below_depth_list_total[i] = [x + y for x, y in zip(prop_below_depth_list_total[i], prop_below_depth_list[i])]
        
        amplicons_length_list_total = [x + y for x, y in zip(amplicons_length_list_total, amplicons_length_list)]
        dropout_lengths_list_total = [x + y for x, y in zip(dropout_lengths_list_total, dropout_lengths_list)]
        
    
    
    storing_list_total = [median_depth_list_total, mean_depth_list_total, het_sites_list_total, prop_below_depth_list_total, length_list_total, amplicons_length_list_total, dropout_lengths_list_total]
    
    # Save the total median depth list
    with open(final_path, "wb") as f:
        pickle.dump(storing_list_total, f)

    print("Done!")
