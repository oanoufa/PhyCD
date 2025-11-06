source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 6_fcc_m -t 24:00:00 --mem=256G \
    -o /nfs/research/goldman/anoufa/data/out_err/6/6_fcc_m.out \
    -e /nfs/research/goldman/anoufa/data/out_err/6/6_fcc_m.err \
    --wrap="${pypy_path} \
        /nfs/research/goldman/anoufa/src/dpca/6_find_contaminants_candidates.py \
        --masked_or_random masked"

sbatch -J 6_fcc_r -t 24:00:00 --mem=256G \
    -o /nfs/research/goldman/anoufa/data/out_err/6/6_fcc_r.out \
    -e /nfs/research/goldman/anoufa/data/out_err/6/6_fcc_r.err \
    --wrap="${pypy_path} \
        /nfs/research/goldman/anoufa/src/dpca/6_find_contaminants_candidates.py \
        --masked_or_random random"