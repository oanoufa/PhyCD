import requests
import shutil
import os
import time

if __name__ == "__main__":

    sample_list = ['SRR17245202']

    for read_name in sample_list:
        print(f"Processing {read_name}")
        url = "https://www.ebi.ac.uk/ena/portal/api/search"
        params = {
            "result": "read_run",
            "query": "run_accession=" + read_name,
            "fields": "fastq_ftp",
            "limit": 1
        }

        response = requests.get(url, params=params)
        lines = response.text.strip().split("\n")
        result = lines[1]
        ftp_urls = result.split("\t")[1]

        ftp_urls = ftp_urls.split(";") 
        target_dir = "/nfs/research/goldman/anoufa/data/samples/" + read_name + "/"
        os.makedirs(target_dir, exist_ok=True)

        for link in ftp_urls:
            filename = os.path.basename(link)
            url = "https://" + link
            out_path = os.path.join(target_dir, filename)
            if os.path.exists(out_path):
                print(f"File {filename} already exists at {out_path}, skipping download.")
                continue
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(out_path, "wb") as f:
                    shutil.copyfileobj(r.raw, f)
            print(f"Downloaded {filename} to {out_path}")
        
        # Wait 2s
        time.sleep(2)
        
    print("All downloads completed.")
    
    # Write the sample list to the txt file
    sample_names_path = "/nfs/research/goldman/anoufa/data/samples/z_logs/sample_names.txt"
    with open(sample_names_path, "w") as f:
        for read_name in sample_list:
            f.write(read_name + "\n")
            
    print("Wrote samples names to txt file")
    print("Sending read mapping jobs to the cluster now...")
    os.system("sbatch /nfs/research/goldman/anoufa/shell_scripts/all_fastq_to_bam.slurm")
    print("All jobs sent.")