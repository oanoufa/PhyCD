source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J fd_eyre -t 02:00:00 --mem=64G \
    -o /nfs/research/goldman/anoufa/data/out_err/fd_eyre.out \
    -e /nfs/research/goldman/anoufa/data/out_err/fd_eyre.err \
    --wrap="${pyenv_path} /nfs/research/goldman/anoufa/src/dpca/format_data_eyre.py"
