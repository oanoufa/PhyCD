#!/bin/bash
wrap_cmd="/nfs/research/goldman/anoufa/.venv/bin/pypy3.10 /nfs/research/goldman/anoufa/src/hmm/find_cont_samples.py --n_batch 128 --batch_id \$i --n_het_thr 3"

time="02:00:00"
mem="4G"

for i in $(seq 0 127); do
  if [ "$i" -lt 5 ]; then
    sbatch -J fcs_$i -t $time --mem=$mem \
      -o /nfs/research/goldman/anoufa/data/dpca/out_err/fcs_$i.out \
      -e /nfs/research/goldman/anoufa/data/dpca/out_err/fcs_$i.err \
      --wrap="$(eval echo "$wrap_cmd")"
  else
    sbatch -J fcs_$i -t $time --mem=$mem \
      -o /dev/null \
      -e /dev/null \
      --wrap="$(eval echo "$wrap_cmd")"
  fi
done