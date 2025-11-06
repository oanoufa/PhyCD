#!/bin/bash
time="02:00:00"
mem="4G"

source /nfs/research/goldman/anoufa/shell_scripts/params.sh
n_batch=1
for i in $(seq 0 $(($n_batch - 1))); do

  wrap_cmd="${pyenv_path} \
    /nfs/research/goldman/anoufa/src/dpca/HMM_playground.py \
    --batch_id ${i}"
  sbatch -J hmm_$i -t $time --mem=$mem \
    -o /nfs/research/goldman/anoufa/data/out_err/hmm_$i.out \
    -e /nfs/research/goldman/anoufa/data/out_err/hmm_$i.err \
    --wrap="$(eval echo "$wrap_cmd")"
done
