source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J play -t 02:00:00 --mem=64G \
    -o /nfs/research/goldman/anoufa/data/out_err/pa.out \
    -e /nfs/research/goldman/anoufa/data/out_err/pa.err \
    --wrap="${pyenv_path} /nfs/research/goldman/anoufa/security_side_project/src/process_access.py"
