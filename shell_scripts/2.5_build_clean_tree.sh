source /nfs/research/goldman/anoufa/shell_scripts/params.sh

input_alignment_new="/nfs/research/goldman/anoufa/data/MAPLE_input/alignment_files/clean_tree_alignment_file_${param_term}.maple"
input_alignment_old="/nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_tree_alignment_file_${param_term}.maple"
job_name="ct_${param_term}"
inputTree=/nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_0.1_3_round2_preliminary_tree.tree

sbatch -J ${job_name} -t 300:00:00 --mem=100G \
    -o /nfs/research/goldman/anoufa/data/out_err/2/${job_name}.out \
    -e /nfs/research/goldman/anoufa/data/out_err/2/${job_name}.err \
    --wrap="${pypy_path} \
        ${MAPLE_path} \
        --rateVariation \
        --model UNREST \
        --input ${input_alignment_old} \
        --output /nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_aug${param_term} \
        --overwrite \
        --inputTree ${inputTree}"
        # --inputRates ${inputRates}"
        # --noFastTopologyInitialSearch \
        # --numCores ${num_cores} \
        # --numTopologyImprovements 1"

# Copy the alignment file to the output directory for record-keeping
# cp ${input_alignment_new} \
#    ${input_alignment_old}