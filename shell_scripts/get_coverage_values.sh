source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J get_c_val -t 02:00:00 --mem=64G \
    -o /nfs/research/goldman/anoufa/data/out_err/get_c_val.out \
    -e /nfs/research/goldman/anoufa/data/out_err/get_c_val.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/get_coverage_values.py"