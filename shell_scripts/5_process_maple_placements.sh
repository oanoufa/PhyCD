source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 5_pMSP_m -t 04:00:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/5/5_pMSP_m.out \
    -e /nfs/research/goldman/anoufa/data/out_err/5/5_pMSP_m.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/5_process_maple_placements.py \
        --masked_or_random masked"

sbatch -J 5_pMSP_r -t 04:00:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/5/5_pMSP_r.out \
    -e /nfs/research/goldman/anoufa/data/out_err/5/5_pMSP_r.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/5_process_maple_placements.py \
        --masked_or_random random"