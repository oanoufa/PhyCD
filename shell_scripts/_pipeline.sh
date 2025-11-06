#!/bin/bash

# Change params before starting this in /nfs/research/goldman/anoufa/src/dpca/params.py

# Then run this file like this:
# nohup bash /nfs/research/goldman/anoufa/shell_scripts/_pipeline.sh > /nfs/research/goldman/anoufa/data/dpca/out_err/pipeline.out 2> /nfs/research/goldman/anoufa/data/dpca/out_err/pipeline.err &

# Load snakemake module if not already loaded
if ! command -v snakemake &> /dev/null; then
  module load snakemake/8.5.2
fi


snakemake --executor slurm \
    --jobs 3
# echo "Starting 1_gen_maple_file..."
# bash /nfs/research/goldman/anoufa/shell_scripts/1_gen_maple_file.sh

# WAIT_TIME_1=7200  # 2 hours
# echo "Script 1 sent. Waiting for $WAIT_TIME_1 seconds..."
# sleep $WAIT_TIME_1

# echo "Starting 2_concat_gmf_output..."
# bash /nfs/research/goldman/anoufa/shell_scripts/2_concat_gmf_output.sh

# WAIT_TIME_2=600 # 10 minutes
# echo "Script 2 sent. Waiting for $WAIT_TIME_2 seconds..."
# sleep $WAIT_TIME_2

# echo "Starting 3_maple_sample_placement..."
# bash /nfs/research/goldman/anoufa/shell_scripts/3_maple_sample_placement.sh

# WAIT_TIME_3=50400 # 14 hours
# echo "Script 3 sent. Waiting for $WAIT_TIME_3 seconds..."
# sleep $WAIT_TIME_3

# echo "Starting 4_concat_maple_output..."
# bash /nfs/research/goldman/anoufa/shell_scripts/4_concat_maple_output.sh

# WAIT_TIME_4=2700 # 45 minutes
# echo "Script 4 sent. Waiting for $WAIT_TIME_4 seconds..."
# sleep $WAIT_TIME_4

# echo "Starting 5_handle_sample_placements..." # This takes about 3 hours
# bash /nfs/research/goldman/anoufa/shell_scripts/5_handle_sample_placements.sh
# echo "Script 5 sent."

# echo "All scripts have been sent."