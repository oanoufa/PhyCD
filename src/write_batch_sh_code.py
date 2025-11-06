import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Detect potential contaminated areas depending on several thresholds applied on the Viridian samples.')

parser.add_argument('--nbatch', type=int, default=256,
                    help='Number of batches to generate. Default is 256.')
parser.add_argument('--data_storing_dir', type=str, default="/nfs/research/goldman/anoufa/data/dpca/",
                    help='Data storing directory.')

parser.add_argument('--het_thr', type=float, default=0.9,
                    help='Threshold for the heterozygous proportion. Default is 0.9.')
parser.add_argument('--depth_thr', type=float, default=0.2,
                    help='Threshold for the depth. Default is 0.2 * median depth.')
parser.add_argument('--prop_under_depth_thr', type=float, default=0.3,
                    help='Threshold for the maximal proportion of positions with clean depth below the threshold. Default is 0.3.')
parser.add_argument('--typical_depth_thr', type=int, default=800,
                    help='Median depth threshold. Default is 800.')
parser.add_argument('--path_ref_seq', type=str, default="/nfs/research/goldman/anoufa/data/NC_045512.2.fasta",
                    help='Path to the reference sequence. Default is /nfs/research/goldman/anoufa/data/NC_045512.2.fasta.')
args = parser.parse_args()

# Get the arguments
n_batch = args.nbatch
data_storing_dir = args.data_storing_dir

het_thr = args.het_thr
depth_thr = args.depth_thr
prop_under_depth_thr = args.prop_under_depth_thr
typical_depth_thr = args.typical_depth_thr
path_ref_seq = args.path_ref_seq

batches_dir = os.path.join(data_storing_dir, "batches")
os.makedirs(data_storing_dir, exist_ok=True)

# Clean up and recreate batches directory
if os.path.exists(batches_dir):
    shutil.rmtree(batches_dir)
os.makedirs(batches_dir)

if __name__ == "__main__":
    
    bash_script_path = os.path.join(data_storing_dir, "dpca_batch.sh")
    # Generate the shell script to run the batches
    print(f"Writing SLURM batch script to: {bash_script_path}")
    with open(bash_script_path, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write(f"for i in $(seq 0 {n_batch - 1}); do\n")
        f.write("  sbatch -J dpca_$i -t 48:00:00 --mem=8G \\\n")
        f.write(f"    -o {data_storing_dir}/out_err/dpca_batch_$i.out \\\n")
        f.write(f"    -e {data_storing_dir}/out_err/dpca_batch_$i.err \\\n")
        f.write("    --wrap=\"/nfs/research/goldman/anoufa/.venv/bin/pypy3.10 /nfs/research/goldman/anoufa/src/dpca.py "
            + f"--n_batch {n_batch} "
            + "--batch_id $i "
            + f"--het_thr {het_thr} "
            + f"--depth_thr {depth_thr} "
            + f"--prop_under_depth_thr {prop_under_depth_thr} "
            + f"--typical_depth_thr {typical_depth_thr} "
            + f"--path_ref_seq {path_ref_seq} "
            + "\"\n")
        f.write("done\n")

    print("Done!")
    