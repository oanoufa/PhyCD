source /nfs/research/goldman/anoufa/shell_scripts/params.sh

sbatch -J 4_cmo -t 01:00:00 --mem=16G \
    -o /nfs/research/goldman/anoufa/data/out_err/4/4_cmo.out \
    -e /nfs/research/goldman/anoufa/data/out_err/4/4_cmo.err \
    --wrap="${pyenv_path} \
        /nfs/research/goldman/anoufa/src/dpca/4_concat_maple_output.py"
