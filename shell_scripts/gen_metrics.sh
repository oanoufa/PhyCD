#!/bin/bash

n_batch=256

for i in $(seq 0 $(($n_batch - 1))); do

  wrap_cmd="/nfs/research/goldman/anoufa/.venv/bin/pypy3.10 /nfs/research/goldman/anoufa/src/dpca/gen_metrics.py --n_batch $n_batch --batch_id ${i} --het_thrs 0.7 0.8 0.9 0.95 --depth_thrs 0.2 0.1 0.05 0.02 0.01"

  if [ "$i" -lt 5 ]; then
    sbatch -J genMet_$i -t 01:30:00 --mem=4G \
      -o /nfs/research/goldman/anoufa/data/dpca/out_err/genMet_$i.out \
      -e /nfs/research/goldman/anoufa/data/dpca/out_err/genMet_$i.err \
      --wrap="$(eval echo "$wrap_cmd")"
  else
    sbatch -J genMet_$i -t 01:30:00 --mem=4G \
      -o /dev/null \
      -e /dev/null \
      --wrap="$(eval echo "$wrap_cmd")"
  fi
done
