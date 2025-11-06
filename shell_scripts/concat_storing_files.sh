sbatch -J CSF -t 01:30:00 --mem=2G \
    -o /nfs/research/goldman/anoufa/data/dpca/out_err/CSF.out \
    -e /nfs/research/goldman/anoufa/data/dpca/out_err/CSF.err \
    --wrap="/nfs/research/goldman/anoufa/.venv/bin/pypy3.10 /nfs/research/goldman/anoufa/src/dpca/concat_storing_files.py"
