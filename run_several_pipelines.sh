#!/bin/bash
#SBATCH --job-name=phycd_sm_array
#SBATCH --output=/nfs/research/goldman/anoufa/pipeline/snakemake_files/snake_master_%a.out
#SBATCH --error=/nfs/research/goldman/anoufa/pipeline/snakemake_files/snake_master_%a.err
#SBATCH --time=720:00:00
#SBATCH --mem=32G
#SBATCH --array=0-0

echo "Loading snakemake module..."

# Load Snakemake if not available
if ! command -v snakemake &> /dev/null; then
  module load snakemake/8.5.2
fi

echo "Snakemake module loaded."

# Parameter grid
depth_thr_list=(0.05 0.1 0.15 0.2)
max_dropout_masked_list=(0 0 0 0)
het_thr_list=(0.1 0.1 0.1 0.1)
max_het_list=(3 3 3 3)

# Select based on job ID
depth_thr=${depth_thr_list[$SLURM_ARRAY_TASK_ID]}
het_thr=${het_thr_list[$SLURM_ARRAY_TASK_ID]}
max_n_het_sites=${max_het_list[$SLURM_ARRAY_TASK_ID]}
max_dropout_masked=${max_dropout_masked_list[$SLURM_ARRAY_TASK_ID]}

path_to_snakefile="/nfs/research/goldman/anoufa/Snakefile"

echo "Unlocking workflow..."
snakemake -s $path_to_snakefile --unlock

echo "Running pipeline $SLURM_ARRAY_TASK_ID with:"
echo "depth_thr=$depth_thr | max_dropout_masked=$max_dropout_masked | het_thr=$het_thr | max_n_het_sites=$max_n_het_sites"

snakemake -s $path_to_snakefile \
          --jobs 512 \
          --cores 14 \
          --quiet rules \
          --latency-wait 30 \
          --config depth_thr=${depth_thr} max_dropout_masked=${max_dropout_masked} het_thr=${het_thr} max_n_het_sites=${max_n_het_sites} \
          --until concat_eyre_output \
          --executor slurm 

