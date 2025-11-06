source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 8_hhmm_m -t 4:00:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/8/8_hhmm_m.out \
    -e /nfs/research/goldman/anoufa/data/out_err/8/8_hhmm_m.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/8_heterozygosity_hmm.py \
        --masked_or_random masked"

sbatch -J 8_hhmm_r -t 4:00:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/8/8_hhmm_r.out \
    -e /nfs/research/goldman/anoufa/data/out_err/8/8_hhmm_r.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/8_heterozygosity_hmm.py \
        --masked_or_random random"