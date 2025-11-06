source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 7_fff -t 00:30:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/7/7_fff.out \
    -e /nfs/research/goldman/anoufa/data/out_err/7/7_fff.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/7_figs.py"