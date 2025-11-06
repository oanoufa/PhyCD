sbatch -J gcs -t 03:30:00 --mem=32G \
    -o /nfs/research/goldman/anoufa/data/dpca/out_err/gcs.out \
    -e /nfs/research/goldman/anoufa/data/dpca/out_err/gcs.err \
    --wrap="/homes/anoufa/.pyenv/versions/3.11.6/bin/python /nfs/research/goldman/anoufa/src/dpca/get_contaminant_sample.py"
