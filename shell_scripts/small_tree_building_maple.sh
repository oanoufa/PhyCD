
sbatch -J MTB_s_t -t 500:00:00 --mem=100G \
    -o /nfs/research/goldman/anoufa/data/dpca/out_err/MTB_s_t.out \
    -e /nfs/research/goldman/anoufa/data/dpca/out_err/MTB_s_t.err \
    --wrap="/hps/software/users/goldman/pypy3/pypy3.10-v7.3.15-linux64/bin/pypy3.10 /nfs/research/goldman/anoufa/src/dpca/MAPLEv0.7.4.7.py --rateVariation --model UNREST --input /nfs/research/goldman/anoufa/data/MAPLE_input/no_masking_alignment_file_1_0.1_1_0_30000.maple --output /nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_tree --overwrite"
