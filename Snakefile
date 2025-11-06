########################################
# Snakefile - Automated DPCA Pipeline
########################################

# Configuration and parameters

# Number of batches
N_BATCHES = 3 
BATCHES = list(range(N_BATCHES))

# Path to the main folder that will contain all data and results
DATA_DIR = "/nfs/research/goldman/anoufa/data/pipeline"

# Path to the file containing the sample names and paths to their Viridian output folders
SAMPLES_DIR = "/nfs/research/zi/mhunt/Viridian_wf_paper/Vdn_all_ena/Final_archiving/run2viridian_dir.tsv.xz"

# Define all output files for final rule
rule all:
    input:
        # Final outputs from step 7
        "/nfs/research/goldman/anoufa/data/pipeline/processed_results.tsv"

# Step 1: Generate MAPLE files (batched)
rule gen_maple_file:
    output:
        "{DATA_DIR}/batches/maple_alignment_batch_{wildcards.batch}.done"
    resources:
        mem_mb=4000,
        time="02:00:00"
    shell:
        r"""
        source /nfs/research/goldman/anoufa/shell_scripts/params.sh
        $pyenv_path /nfs/research/goldman/anoufa/src/dpca/1_gen_maple_file.py \
            --batch_id {wildcards.batch} \
            --data_dir {DATA_DIR} \
            --samples_dir {SAMPLES_DIR} \
            > data/out_err/1/1_{wildcards.batch}.out 2> data/out_err/1/1_{wildcards.batch}.err
        """

# Step 2: Process GMF output (after all batches complete)
rule process_gmf_output:
    input:
        expand("{DATA_DIR}/1_batches/maple_alignment_batch_{batch}.done", batch=BATCHES)
    output:
        "{DATA_DIR}/processed_results.tsv"
    resources:
        mem_mb=32000,
        time="02:00:00"
    shell:
        r"""
        source /nfs/research/goldman/anoufa/shell_scripts/params.sh
        $pyenv_path /nfs/research/goldman/anoufa/src/dpca/2_process_gmf_output.py \
            > data/out_err/2/2_cmf.out 2> data/out_err/2/2_cmf.err
        """

# Step 3: MAPLE sample placement (batched)
rule maple_sample_placement:
    resources:
        mem_mb=48000,
        time="09:00:00"
    shell:
        r"""
        source /nfs/research/goldman/anoufa/shell_scripts/params.sh
        
        input_file="data/MAPLE_input/1_batches/maple_alignment_${{masking_method}}_batch{wildcards.batch}_${{param_term}}.maple"
        output_dir="data/MAPLE_output/batches/output_{wildcards.batch}"
        
        $pypy_path $MAPLE_path \
            --inputTree $inputTree \
            --input $input_file \
            --findSamplePlacement \
            --model=UNREST \
            --rateVariation \
            --minBranchSupport 0.005 \
            --inputRates $inputRates \
            --output $output_dir \
            --overwrite \
            > data/out_err/3/3_{wildcards.batch}.out 2> data/out_err/3/3_{wildcards.batch}.err
        """

# # Step 4: Concatenate MAPLE output
# rule concat_maple_output:
#     input:
#         expand("data/MAPLE_output/batches/output_{batch}/sampleLikelihoodPlacement.txt", batch=BATCHES)
#     output:
#         "data/out_err/4/4_cmo.out",
#         "results/maple_concat.done"
#     resources:
#         mem_mb=16000,
#         time="01:00:00"
#     shell:
#         r"""
#         source /nfs/research/goldman/anoufa/shell_scripts/params.sh
#         mkdir -p data/out_err/4
#         $pyenv_path /nfs/research/goldman/anoufa/src/dpca/4_concat_maple_output.py \
#             > {output[0]} 2>&1
#         touch {output[1]}
#         """

# # Step 5: Process MAPLE placements (masked variant)
# rule process_maple_placements_masked:
#     input:
#         "results/maple_concat.done"
#     output:
#         "data/out_err/5/5_pMSP_m.out",
#         "results/placements_masked.done"
#     resources:
#         mem_mb=32000,
#         time="04:00:00"
#     shell:
#         r"""
#         source /nfs/research/goldman/anoufa/shell_scripts/params.sh
#         mkdir -p data/out_err/5
#         $pyenv_path /nfs/research/goldman/anoufa/src/dpca/5_process_maple_placements.py \
#             --masked_or_random masked \
#             > {output[0]} 2>&1
#         touch {output[1]}
#         """

# # Step 5: Process MAPLE placements (random variant)
# rule process_maple_placements_random:
#     input:
#         "results/maple_concat.done"
#     output:
#         "data/out_err/5/5_pMSP_r.out",
#         "results/placements_random.done"
#     resources:
#         mem_mb=32000,
#         time="04:00:00"
#     shell:
#         r"""
#         source /nfs/research/goldman/anoufa/shell_scripts/params.sh
#         mkdir -p data/out_err/5
#         $pyenv_path /nfs/research/goldman/anoufa/src/dpca/5_process_maple_placements.py \
#             --masked_or_random random \
#             > {output[0]} 2>&1
#         touch {output[1]}
#         """

# # Step 6: Find contaminant candidates (masked variant)
# rule find_contaminant_candidates_masked:
#     input:
#         "results/placements_masked.done"
#     output:
#         "data/out_err/6/6_fcc_m.out",
#         "results/contaminants_masked.done"
#     resources:
#         mem_mb=256000,
#         time="24:00:00"
#     shell:
#         r"""
#         source /nfs/research/goldman/anoufa/shell_scripts/params.sh
#         mkdir -p data/out_err/6
#         $pypy_path /nfs/research/goldman/anoufa/src/dpca/6_find_contaminants_candidates.py \
#             --masked_or_random masked \
#             > {output[0]} 2>&1
#         touch {output[1]}
#         """

# # Step 6: Find contaminant candidates (random variant)
# rule find_contaminant_candidates_random:
#     input:
#         "results/placements_random.done"
#     output:
#         "data/out_err/6/6_fcc_r.out",
#         "results/contaminants_random.done"
#     resources:
#         mem_mb=256000,
#         time="24:00:00"
#     shell:
#         r"""
#         source /nfs/research/goldman/anoufa/shell_scripts/params.sh
#         mkdir -p data/out_err/6
#         $pypy_path /nfs/research/goldman/anoufa/src/dpca/6_find_contaminants_candidates.py \
#             --masked_or_random random \
#             > {output[0]} 2>&1
#         touch {output[1]}
#         """

# # Step 7: Generate figures (final step)
# rule generate_figures:
#     input:
#         "results/contaminants_masked.done",
#         "results/contaminants_random.done"
#     output:
#         "data/out_err/7/7_fff.out",
#         "figs/final_results.done"
#     resources:
#         mem_mb=32000,
#         time="00:30:00"
#     shell:
#         r"""
#         source /nfs/research/goldman/anoufa/shell_scripts/params.sh
#         mkdir -p data/out_err/7 figs
#         $pyenv_path /nfs/research/goldman/anoufa/src/dpca/7_figs.py \
#             > {output[0]} 2>&1
#         touch {output[1]}
#         """
