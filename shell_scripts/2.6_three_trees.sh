
source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 2_6_mas -t 720:00:00 --mem=200G \
    -o /nfs/research/goldman/anoufa/data/out_err/2/2_6_mas.out \
    -e /nfs/research/goldman/anoufa/data/out_err/2/2_6_mas.err \
    --wrap="${pypy_path} ${MAPLE_path} \
    --rateVariation --model UNREST \
    --input /nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/masked_alignment_${param_term}.maple \
    --output /nfs/research/goldman/anoufa/data/MAPLE_output/masked/mas \
    --overwrite --noFastTopologyInitialSearch --numTopologyImprovements 0"


sbatch -J 2_6_unma -t 720:00:00 --mem=200G \
    -o /nfs/research/goldman/anoufa/data/out_err/2/2_6_unma.out \
    -e /nfs/research/goldman/anoufa/data/out_err/2/2_6_unma.err \
    --wrap="${pypy_path} ${MAPLE_path} \
    --rateVariation --model UNREST \
    --input /nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/unmasked_alignment_${param_term}.maple \
    --output /nfs/research/goldman/anoufa/data/MAPLE_output/unmasked/unma \
    --overwrite --noFastTopologyInitialSearch --numTopologyImprovements 0"

    
sbatch -J 2_6_rand -t 720:00:00 --mem=200G \
    -o /nfs/research/goldman/anoufa/data/out_err/2/2_6_rand.out \
    -e /nfs/research/goldman/anoufa/data/out_err/2/2_6_rand.err \
    --wrap="${pypy_path} ${MAPLE_path} \
    --rateVariation --model UNREST \
    --input /nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/random_alignment_${param_term}.maple \
    --output /nfs/research/goldman/anoufa/data/MAPLE_output/randomly_masked/rand \
    --overwrite --noFastTopologyInitialSearch --numTopologyImprovements 0"