########################################
# Snakefile pipeline
########################################

# Configuration and parameters

# Number of batches
N_BATCHES = 512
BATCHES = list(range(N_BATCHES))

# Path to the main folder that will contain all data and results
DATA_DIR = "/nfs/research/goldman/anoufa/pipeline"

# Path to the python scripts
SCRIPTS_DIR = "/nfs/research/goldman/anoufa/src/dpca"

# Logs and output directory
OUT_ERR_DIR = "/nfs/research/goldman/anoufa/data/out_err"

# Path to the file containing the sample names and paths to their Viridian output folders
SAMPLES_DIR = "/nfs/research/goldman/anoufa/pipeline/run2viridian_dir.tsv.xz"

# Path to the metadata file containing all sample metadata
METADATA = "/nfs/research/goldman/anoufa/data/MAPLE_input/others/viridian_samples.metadata.tsv"


# Path to the parameters file where paths and parameters will be saved during the pipeline
PARAM_PY_PATH = "/nfs/research/goldman/anoufa/src/dpca/_params.py"
PARAM_SH_PATH = PARAM_PY_PATH.replace(".py", ".sh")

# Path to the reference sequence used for the pipeline
path_ref_seq = "/nfs/research/goldman/anoufa/data/MAPLE_input/maple_ref_lower.fasta"


# Paths to interpreters and MAPLE script
PYENV_PATH="/homes/anoufa/.pyenv/versions/3.11.6/bin/python3"
PYPY_PATH="/nfs/research/goldman/anoufa/.venv/bin/pypy3.10"
MAPLE_PATH="/nfs/research/goldman/anoufa/src/dpca/MAPLEv0.7.5.py"

# Parameters of the pipeline
# MASKING THR, works ~the same with 0.15 and 0.2 (positions are masked if their coverage is under median depth * masking_thr)
depth_thr = 0.1
# HET_THR, positions are considered heterozygous if the proportion of the consensus allele is under 1 - het_thr (= max prop of minor alleles)
het_thr = 0.1
# MAX NUMBER OF HET SITES ALLOWED TO CONSIDER A SAMPLE CLEAN ENOUGH TO BE PART OF THE CLEAN TREE (and excluded from the analysis)
max_n_het_sites = 3
param_term = f"{depth_thr}_{max_n_het_sites}"

# Define all output files for final rule
rule all:
    input:
        # Final outputs from step 7
        f"{DATA_DIR}/done_files/7_figs.done"

# Step 1: Generate MAPLE files (batched)
rule gen_maple_file:
    output:
        f"{DATA_DIR}/done_files/1_maple_alignment_batch_{{batch}}.done"
    resources:
        mem_mb=4000,
        runtime=120 # in minutes
    shell:
        """
        mkdir -p {DATA_DIR}/done_files
        mkdir -p {DATA_DIR}/1
        mkdir -p {OUT_ERR_DIR}/1
        {PYENV_PATH} {SCRIPTS_DIR}/1_gen_maple_file.py \
            --n_batch {N_BATCHES} \
            --batch_id {wildcards.batch} \
            --data_dir {DATA_DIR} \
            --samples_dir {SAMPLES_DIR} \
            --param_path {PARAM_PY_PATH} \
            --depth_thr {depth_thr} \
            --max_n_het_sites {max_n_het_sites} \
            --path_ref_seq {path_ref_seq} \
            --het_thr {het_thr} \
            > {OUT_ERR_DIR}/1/1_{wildcards.batch}.out 2> {OUT_ERR_DIR}/1/1_{wildcards.batch}.err
        """

# Step 2: Process GMF output (after all batches complete)
rule process_gmf_output:
    input:
        expand(f"{DATA_DIR}/done_files/1_maple_alignment_batch_{{batch}}.done", batch=BATCHES)
    output:
        f"{DATA_DIR}/done_files/2_process_gmf_output.done"
    resources:
        mem_mb=32000,
        runtime=120 # in minutes
    shell:
        """
        mkdir -p {DATA_DIR}/2
        mkdir -p {OUT_ERR_DIR}/2
        source {PARAM_SH_PATH}
        {PYENV_PATH} {SCRIPTS_DIR}/2_process_gmf_output.py \
            > {OUT_ERR_DIR}/2/2_cmf.out 2> {OUT_ERR_DIR}/2/2_cmf.err
        """

# Step 3a: Generate clean tree
rule generate_clean_tree:
    input:
        f"{DATA_DIR}/done_files/2_process_gmf_output.done"
    output:
        f"{DATA_DIR}/done_files/3_generate_clean_tree.done"
    resources:
        mem_mb=64000,
        runtime=12000 # in minutes
    shell:
        """
        mkdir -p {DATA_DIR}/3/clean_tree
        mkdir -p {OUT_ERR_DIR}/3
        source {PARAM_SH_PATH}

        {PYPY_PATH} {MAPLE_PATH} \
        --rateVariation \
        --model UNREST \
        --input ${{final_clean_tree_path}} \
        --output {DATA_DIR}/3/clean_tree/clean_${{param_term}} \
        --overwrite \
        > {OUT_ERR_DIR}/3/3_gct.out 2> {OUT_ERR_DIR}/3/3_gct.err

        touch {output}
        """

# Step 3b: Generate JSONL tree for Taxonium visualization
rule generate_taxonium_tree:
    input:
        f"{DATA_DIR}/done_files/3_generate_clean_tree.done"
    output:
        f"{DATA_DIR}/done_files/3_newick_to_taxonium.done"
    resources:
        mem_mb=32000,
        runtime=240  # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        inputTree="{DATA_DIR}/3/clean_tree/clean_${{param_term}}_tree.tree"

        {PYENV_PATH} -m taxoniumtools.newick_to_taxonium \
            -i ${{inputTree}} \
            -m {METADATA} \
            -o {DATA_DIR}/3/clean_tree/clean_${{param_term}}_tree.jsonl \
            -c Country,Viridian_pangolin_1.29 \
            --key_column Run \
            > {OUT_ERR_DIR}/3/3_ntt.out 2> {OUT_ERR_DIR}/3/3_ntt.err

        touch {output}
        """

# Step 3c: MAPLE sample placement (batched)
rule maple_sample_placement:
    input:
        f"{DATA_DIR}/done_files/3_newick_to_taxonium.done"
    output:
        f"{DATA_DIR}/done_files/3_maple_sample_placement_{{batch}}.done"
    resources:
        mem_mb=48000,
        runtime=540 # in minutes
    shell:
        """
        mkdir -p {DATA_DIR}/3
        mkdir -p {OUT_ERR_DIR}/3
        source {PARAM_SH_PATH}


        inputTree="{DATA_DIR}/3/clean_tree/clean_${{param_term}}_tree.tree"
        inputRates="{DATA_DIR}/3/clean_tree/clean_${{param_term}}_subs.txt"
        inputFile="{DATA_DIR}/1/maple_alignment_batch{wildcards.batch}_${{param_term}}.maple"

        {PYPY_PATH} {MAPLE_PATH} \
            --inputTree ${{inputTree}} \
            --input ${{inputFile}} \
            --findSamplePlacement \
            --model UNREST \
            --rateVariation \
            --minBranchSupport 0.005 \
            --inputRates ${{inputRates}} \
            --output "{DATA_DIR}/3/output_{wildcards.batch}" \
            --overwrite \
            > {OUT_ERR_DIR}/3/3_{wildcards.batch}.out 2> {OUT_ERR_DIR}/3/3_{wildcards.batch}.err

        touch {output}
        """

# Step 4: Concatenate MAPLE output
rule concat_maple_output:
    input:
        expand(f"{DATA_DIR}/done_files/3_maple_sample_placement_{{batch}}.done", batch=BATCHES)
    output:
        f"{DATA_DIR}/done_files/4_concat_maple_output.done"
    resources:
        mem_mb=16000,
        runtime=60 # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        mkdir -p {DATA_DIR}/4
        mkdir -p {OUT_ERR_DIR}/4
        {PYENV_PATH} {SCRIPTS_DIR}/4_concat_maple_output.py \
            > {OUT_ERR_DIR}/4/4_cmo.out 2> {OUT_ERR_DIR}/4/4_cmo.err
        """

# Step 5a: Process MAPLE placements (masked variant)
rule process_maple_placements_masked:
    input:
        f"{DATA_DIR}/done_files/4_concat_maple_output.done"
    output:
        f"{DATA_DIR}/done_files/5_process_maple_placements_masked.done"
    resources:
        mem_mb=32000,
        runtime=240 # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        mkdir -p {DATA_DIR}/5
        mkdir -p {OUT_ERR_DIR}/5
        {PYENV_PATH} {SCRIPTS_DIR}/5_process_maple_placements.py \
            --masked_or_random masked \
            > {OUT_ERR_DIR}/5/5_pMSP_m.out 2> {OUT_ERR_DIR}/5/5_pMSP_m.err
        """

# Step 5b: Process MAPLE placements (random variant)
rule process_maple_placements_random:
    input:
        f"{DATA_DIR}/done_files/4_concat_maple_output.done"
    output:
        f"{DATA_DIR}/done_files/5_process_maple_placements_random.done"
    resources:
        mem_mb=32000,
        runtime=240 # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        mkdir -p {DATA_DIR}/5
        mkdir -p {OUT_ERR_DIR}/5
        {PYENV_PATH} {SCRIPTS_DIR}/5_process_maple_placements.py \
            --masked_or_random random \
            > {OUT_ERR_DIR}/5/5_pMSP_r.out 2> {OUT_ERR_DIR}/5/5_pMSP_r.err
        """

# Step 6a: Find contaminant candidates (masked variant)
rule find_contaminant_candidates_masked:
    input:
        f"{DATA_DIR}/done_files/5_process_maple_placements_masked.done"
    output:
        f"{DATA_DIR}/done_files/6_find_contaminants_candidates_masked.done"
    resources:
        mem_mb=128000,
        runtime=800 # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        mkdir -p {OUT_ERR_DIR}/6
        mkdir -p {DATA_DIR}/6
        {PYENV_PATH} {SCRIPTS_DIR}/6_find_contaminants_candidates.py \
            --masked_or_random masked \
            > {OUT_ERR_DIR}/6/6_fcc_m.out 2> {OUT_ERR_DIR}/6/6_fcc_m.err
        """

# Step 6b: Find contaminant candidates (random variant)
rule find_contaminant_candidates_random:
    input:
        f"{DATA_DIR}/done_files/5_process_maple_placements_random.done"
    output:
        f"{DATA_DIR}/done_files/6_find_contaminants_candidates_random.done"
    resources:
        mem_mb=128000,
        runtime=800 # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        mkdir -p {OUT_ERR_DIR}/6
        mkdir -p {DATA_DIR}/6
        {PYENV_PATH} {SCRIPTS_DIR}/6_find_contaminants_candidates.py \
            --masked_or_random random \
            > {OUT_ERR_DIR}/6/6_fcc_r.out 2> {OUT_ERR_DIR}/6/6_fcc_r.err
        """

# Step 7: Generate figures (final step)
rule generate_figures:
    input:
        f"{DATA_DIR}/done_files/6_find_contaminants_candidates_masked.done",
        f"{DATA_DIR}/done_files/6_find_contaminants_candidates_random.done"
    output:
        f"{DATA_DIR}/done_files/7_figs.done"
    resources:
        mem_mb=32000,
        runtime=30 # in minutes
    shell:
        """
        source {PARAM_SH_PATH}
        mkdir -p {OUT_ERR_DIR}/7
        mkdir -p {DATA_DIR}/7
        {PYENV_PATH} {SCRIPTS_DIR}/7_figs.py \
            > {OUT_ERR_DIR}/7/7_figs.out 2> {OUT_ERR_DIR}/7/7_figs.err
        """
