source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 0_sfhmm -t 03:00:00 --mem=64G \
    -o /nfs/research/goldman/anoufa/data/out_err/0/0_sfhmm.out \
    -e /nfs/research/goldman/anoufa/data/out_err/0/0_sfhmm.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/0_setup_fit_hmm.py"

