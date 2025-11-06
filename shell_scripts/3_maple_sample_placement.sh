

time="09:00:00"
mem="48G"

source /nfs/research/goldman/anoufa/shell_scripts/params.sh

for i in $(seq 0 $(($n_batch - 1))); do
  input_file="/nfs/research/goldman/anoufa/data/MAPLE_input/1_batches/maple_alignment_${masking_method}_batch${i}_${param_term}.maple"
  output_dir="/nfs/research/goldman/anoufa/data/MAPLE_output/batches/output_${i}"
  wrap_cmd="${pypy_path} \
    ${MAPLE_path} \
    --inputTree ${inputTree} \
    --input ${input_file} \
    --findSamplePlacement \
    --model=UNREST \
    --rateVariation \
    --minBranchSupport 0.005 \
    --inputRates ${inputRates} \
    --output ${output_dir} \
    --overwrite"


  if [ "$i" -lt 600 ]; then
    sbatch -J 3_$i -t $time --mem=$mem \
      -o /nfs/research/goldman/anoufa/data/out_err/3/3_$i.out \
      -e /nfs/research/goldman/anoufa/data/out_err/3/3_$i.err \
      --wrap="$wrap_cmd"

  # If it's the last job, retain the job ID for dependency
  elif [ "$i" -eq $((n_batch - 1)) ]; then
    jobID=$(sbatch -J 3_$i -t $time --mem=$mem \
      -o /dev/null \
      -e /dev/null \
      --wrap="$wrap_cmd" | awk '{print $4}')

  else
    sbatch -J 3_$i -t $time --mem=$mem \
      -o /dev/null \
      -e /dev/null \
      --wrap="$wrap_cmd"
  fi
done

echo "Submitted all sample placement jobs. Last job ID: $jobID" 

