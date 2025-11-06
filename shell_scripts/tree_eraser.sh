sbatch -J tree_eraser -t 00:05:00 --mem=1G \
  -o /dev/null \
  -e /dev/null \
  --wrap="/nfs/research/goldman/anoufa/.venv/bin/pypy3.10 /nfs/research/goldman/anoufa/src/dpca/tree_eraser.py" 