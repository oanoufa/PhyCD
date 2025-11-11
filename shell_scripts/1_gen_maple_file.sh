#!/bin/bash

time="02:00:00"
mem="4G"

source /nfs/research/goldman/anoufa/shell_scripts/params.sh

for i in $(seq 0 $(($n_batch - 1))); do

  wrap_cmd="${pyenv_path} \
    /nfs/research/goldman/anoufa/src/dpca/1_gen_maple_file.py \
    --batch_id ${i} \
    --data_dir /nfs/research/goldman/anoufa/data/pipeline \
    --samples_dir /nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz \
    --param_path /nfs/research/goldman/anoufa/src/dpca/_params.py"

  sbatch -J 1_$i -t $time --mem=$mem \
    -o /nfs/research/goldman/anoufa/data/out_err/1/1_$i.out \
    -e /nfs/research/goldman/anoufa/data/out_err/1/1_$i.err \
    --wrap="$(eval echo "$wrap_cmd")"
done

# for i in $(seq 0 $(($n_batch - 1))); do

#   wrap_cmd="${pypy_path} \
#     /nfs/research/goldman/anoufa/src/dpca/1_gen_maple_file.py \
#     --batch_id ${i}"

#   if [ "$i" -lt 5 ]; then
#     sbatch -J 1_$i -t $time --mem=$mem \
#       -o /nfs/research/goldman/anoufa/data/dpca/out_err/1_$i.out \
#       -e /nfs/research/goldman/anoufa/data/dpca/out_err/1_$i.err \
#       --wrap="$(eval echo "$wrap_cmd")"
#   else
#     sbatch -J 1_$i -t $time --mem=$mem \
#       -o /dev/null \
#       -e /dev/null \
#       --wrap="$(eval echo "$wrap_cmd")"
#   fi
# done