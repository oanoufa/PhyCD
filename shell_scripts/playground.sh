source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J play -t 02:00:00 --mem=64G \
    -o /nfs/research/goldman/anoufa/data/out_err/play.out \
    -e /nfs/research/goldman/anoufa/data/out_err/play.err \
    --wrap="${pyenv_path} /nfs/research/goldman/anoufa/src/dpca/playground.py"
