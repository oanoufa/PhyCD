#!/bin/bash

# Default parameters
depth_thr=0.05
max_dropout_masked=0
het_thr=0.1
max_n_het_sites=3

path_to_snakefile="./Snakefile"

echo "Unlocking workflow..."
snakemake -s $path_to_snakefile --unlock

echo "Running example pipeline with:"
echo "depth_thr=$depth_thr | max_dropout_masked=$max_dropout_masked | het_thr=$het_thr | max_n_het_sites=$max_n_het_sites"

snakemake -s $path_to_snakefile \
          --cores 1 \
          --latency-wait 5 \
          --config depth_thr=${depth_thr} \
                   max_dropout_masked=${max_dropout_masked} \
                   het_thr=${het_thr} \
                   max_n_het_sites=${max_n_het_sites} \
          --until concat_eyre_output