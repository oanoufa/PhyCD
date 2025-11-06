sbatch -J FMA -t 01:30:00 --mem=16G \
    -o /nfs/research/goldman/anoufa/data/dpca/out_err/FMA.out \
    -e /nfs/research/goldman/anoufa/data/dpca/out_err/FMA.err \
    --wrap="/nfs/research/goldman/anoufa/.venv/bin/pypy3.10 /nfs/research/goldman/anoufa/src/dpca/remove_potentially_cont_from_2M.py"
