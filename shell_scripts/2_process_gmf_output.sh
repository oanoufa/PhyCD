source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 2_cmf -t 02:00:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/out_err/2/2_cmf.out \
    -e /nfs/research/goldman/anoufa/data/out_err/2/2_cmf.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/2_process_gmf_output.py"

