#!/bin/bash

# Change params before starting this in /nfs/research/goldman/anoufa/src/dpca/params.py

# Then run this file like this:
# nohup bash /nfs/research/goldman/anoufa/shell_scripts/_pipeline.sh > /nfs/research/goldman/anoufa/data/dpca/out_err/pipeline.out 2> /nfs/research/goldman/anoufa/data/dpca/out_err/pipeline.err &

# Load snakemake module if not already loaded
if ! command -v snakemake &> /dev/null; then
  module load snakemake/8.5.2
fi

path_to_snakefile="/nfs/research/goldman/anoufa/Snakefile"
n_jobs=512

snakemake -s $path_to_snakefile \
          --jobs $n_jobs \
          --executor slurm \
          --latency-wait 60
          > /nfs/research/goldman/anoufa/data/out_err/snake.out \
          2> /nfs/research/goldman/anoufa/data/out_err/snake.err
