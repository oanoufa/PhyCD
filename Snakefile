########################################
# PhyCD Snakemake pipeline
########################################

# Configuration and parameters
# Number of batches (= number of SLURM jobs for the batched scripts)
n_batch = 512
BATCHES = list(range(n_batch))

# Number of cores to use when multiprocessing
num_cores = 8

# MASKING THR, works ~the same with 0.15 and 0.2 (positions are masked if their coverage is under or equal to median depth * masking_thr) (<=)
depth_thr = float(config.get("depth_thr", 0.1) or 0.1)        # default 0.1 if missing, pM in the paper
# MAX NUMBER OF DROPOUT MASKED POSITIONS ALLOWED TO CONSIDER A SAMPLE CLEAN ENOUGH TO BE PART OF THE CLEAN TREE (and excluded from the analysis) (sample can have max_dropout_masked or less masked sites (<=))
max_dropout_masked = int(config.get("max_dropout_masked", 0) or 0)  # default 0, kappa in the paper
# HET_THR, positions are considered heterozygous if the proportion of the minor allele is over or equal to het_thr (>=)
het_thr = float(config.get("het_thr", 0.1) or 0.1)           # default 0.1, eta in the paper
# MAX NUMBER OF HET SITES ALLOWED TO CONSIDER A SAMPLE CLEAN ENOUGH TO BE PART OF THE CLEAN TREE (and excluded from the analysis), (sample can have max_n_het_sites or less heterozygous sites (<=))
max_n_het_sites = int(config.get("max_n_het_sites", 0) or 0)  # default 0, theta in the paper

# PARAM_TERM is used to keep track of the current parameters and avoid overwriting if parameters are changed
param_term = f"{int(depth_thr*100)}_{max_dropout_masked}_{int(het_thr*100)}_{max_n_het_sites}"
print("Parameters received:", flush=True)
print("depth_thr:", depth_thr, flush=True)
print("max_dropout_masked:", max_dropout_masked, flush=True)
print("het_thr:", het_thr, flush=True)
print("max_n_het_sites:", max_n_het_sites, flush=True)

# 5_process_sample_placements default parameters
# Minimal distance difference between the unmasked and masked placements.
n_diff_mut = 2
# Maximal distance to the tree of the masked placement
masked_max_dist = 5
# Minimal ratio between the proportion of distance reduced and the proportion of genome masked.
masking_ratio = 0

# Path to root dir
root_dir = "/nfs/research/goldman/anoufa"
# Path to the main folder that will contain all data and results
data_dir = f"{root_dir}/pipeline/{param_term}/data"
# Path to the python scripts
scripts_dir = f"{root_dir}/scripts"
# Logs and output directory
out_err_dir = f"{root_dir}/pipeline/{param_term}/out_err"
# Path to the file containing the sample names and paths to their Viridian output folders
samples_dir = f"{root_dir}/input_data/run2viridian_dir.tsv.xz"
# Path to the metadata file containing all sample metadata
metadata = f"{root_dir}/input_data/viridian_samples.metadata.tsv"
# Path to the reference sequence used for the pipeline
path_ref_seq = f"{root_dir}/input_data/maple_ref_lower.fasta"
# Paths to interpreters and MAPLE script
pyenv_path="/homes/anoufa/.pyenv/versions/3.11.6/bin/python3"
pypy_path=f"{root_dir}/.venv/bin/pypy3.10"
maple_path=f"{scripts_dir}/MAPLEv0.7.5.py"
# Compress files or not
compress = 1 # 1 or 0


# Define all output files for final rule
rule all:
    input:
        # Outputs of all steps
        expand(f"{data_dir}/done_files/1/1_maple_alignment_batch_{{batch}}.done", batch=BATCHES),
        f"{data_dir}/done_files/2_process_gmf_output.done",
        f"{data_dir}/done_files/2_generate_clean_tree.done",
        f"{data_dir}/done_files/2_newick_to_taxonium.done",
        expand(f"{data_dir}/done_files/3/3_maple_sample_placement_{{batch}}.done", batch=BATCHES),
        f"{data_dir}/done_files/4_concat_maple_output.done",
        f"{data_dir}/done_files/5_process_maple_placements_masked.done",
        f"{data_dir}/done_files/5_process_maple_placements_random.done",
        f"{data_dir}/done_files/6_find_contaminants_candidates_masked.done",
        f"{data_dir}/done_files/6_find_contaminants_candidates_random.done",
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_random.done",
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_masked.done",
        f"{data_dir}/done_files/8_concat_eyre_output.done"

# Step 1: Generate MAPLE files (batched)
rule gen_maple_file:
    threads: 1
    output:
        f"{data_dir}/done_files/1/1_maple_alignment_batch_{{batch}}.done"
    resources:
        mem_mb=4000,
        runtime=240 # in minutes
    shell:
        """
        mkdir -p {data_dir}/done_files/1
        mkdir -p {data_dir}/1
        mkdir -p {out_err_dir}/1
        {pyenv_path} {scripts_dir}/1_gen_maple_file.py \
            --n_batch {n_batch} \
            --batch_id {wildcards.batch} \
            --data_dir {data_dir} \
            --samples_dir {samples_dir} \
            --depth_thr {depth_thr} \
            --max_n_het_sites {max_n_het_sites} \
            --max_dropout_masked {max_dropout_masked} \
            --path_ref_seq {path_ref_seq} \
            --het_thr {het_thr} \
            --param_term {param_term} \
            --compress {compress} \
            > {out_err_dir}/1/1_{wildcards.batch}_sm.out 2> {out_err_dir}/1/1_{wildcards.batch}_sm.err

        touch {output}
        """

# Step 2a: Process GMF output (after all batches complete)
rule process_gmf_output:
    threads: num_cores
    input:
        expand(f"{data_dir}/done_files/1/1_maple_alignment_batch_{{batch}}.done", batch=BATCHES)
    output:
        f"{data_dir}/done_files/2_process_gmf_output.done"
    resources:
        mem_mb=32000,
        runtime=120 # in minutes
    shell:
        """
        mkdir -p {data_dir}/2
        mkdir -p {out_err_dir}/2
        {pyenv_path} {scripts_dir}/2_process_gmf_output.py \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --path_ref_seq {path_ref_seq} \
            --path_metadata_tsv {metadata} \
            --num_cores {num_cores} \
            --compress {compress} \
            > {out_err_dir}/2/2_cmf_sm.out 2> {out_err_dir}/2/2_cmf_sm.err
        
        touch {output}
        """

# Step 2b: Generate clean tree
rule generate_clean_tree:
    threads: num_cores
    input:
        f"{data_dir}/done_files/2_process_gmf_output.done"
    output:
        f"{data_dir}/done_files/2_generate_clean_tree.done"
    resources:
        mem_mb=1024000,
        runtime=30000 # in minutes = 400 hours = ~21 days
    shell:
        """
        # Add this line if program crashed and you need to start from a checkpoint tree.
        # --inputTree path/to/input/tree.tree \

        # Check if the output file exists
        if [ ! -f "{data_dir}/clean_tree/clean_{param_term}_tree.tree" ]; then
            
            
            mkdir -p {data_dir}/clean_tree

            {pypy_path} {maple_path} \
            --rateVariation \
            --model UNREST \
            --saveInitialTreeEvery 200000 \
            --input {data_dir}/2/alignment_files/clean_tree_alignment_file_{param_term}.maple \
            --output {data_dir}/clean_tree/clean_{param_term} \
            --numCores 14 \
            --overwrite \
            --numTopologyImprovements 1 \
            > {out_err_dir}/2/2_gct_sm.out 2> {out_err_dir}/2/2_gct_sm.err

            touch {output} 
                        
        else
            echo "Clean tree already exists. Skipping rule execution."
            touch {output}
        fi
        """

# Step 2c: Generate JSONL tree for Taxonium visualization
rule generate_taxonium_tree:
    threads: 1
    input:
        f"{data_dir}/done_files/2_generate_clean_tree.done"
    output:
        f"{data_dir}/done_files/2_newick_to_taxonium.done"
    resources:
        mem_mb=32000,
        runtime=240  # in minutes
    shell:
        """
        inputTree="{data_dir}/clean_tree/clean_{param_term}_tree.tree"
        outputTree="{data_dir}/clean_tree/clean_{param_term}_tree.jsonl"

        if [ ! -f ${{outputTree}} ]; then

            {pyenv_path} -m taxoniumtools.newick_to_taxonium \
                -i ${{inputTree}} \
                -m {metadata} \
                -o ${{outputTree}} \
                -c Country,Viridian_pangolin_1.29 \
                --key_column Run \
                > {out_err_dir}/2/2_ntt_sm.out 2> {out_err_dir}/2/2_ntt_sm.err

            touch {output}
        else
            echo "JSONL tree already exists. Skipping rule execution."
            touch {output}
        fi
        """

# Step 3c: MAPLE sample placement (batched)
rule maple_sample_placement:
    threads: 1
    input:
        f"{data_dir}/done_files/2_newick_to_taxonium.done"
    output:
        f"{data_dir}/done_files/3/3_maple_sample_placement_{{batch}}.done"
    resources:
        mem_mb=20000,
        runtime=300 # in minutes = 10 hours
    shell:
        """
        mkdir -p {data_dir}/3
        mkdir -p {out_err_dir}/3
        mkdir -p {data_dir}/done_files/3

        inputTree="{data_dir}/clean_tree/clean_{param_term}_tree.tree"
        inputRates="{data_dir}/clean_tree/clean_{param_term}_subs.txt"
        inputFile="{data_dir}/1/maple_alignment_batch{wildcards.batch}_{param_term}.maple"
        # Append .zst if compress=1
        [ {compress} -eq 1 ] && inputFile="${{inputFile}}.zst"

        {pypy_path} {maple_path} \
            --inputTree ${{inputTree}} \
            --input ${{inputFile}} \
            --findSamplePlacement \
            --model UNREST \
            --numCores 1 \
            --rateVariation \
            --minBranchSupport 0.005 \
            --inputRates ${{inputRates}} \
            --output "{data_dir}/3/output_{wildcards.batch}" \
            --overwrite \
            > {out_err_dir}/3/3_{wildcards.batch}_sm.out 2> {out_err_dir}/3/3_{wildcards.batch}_sm.err

        touch {output}
        """

# Step 4: Concatenate MAPLE output
rule concat_maple_output:
    threads: 1
    input:
        expand(f"{data_dir}/done_files/3/3_maple_sample_placement_{{batch}}.done", batch=BATCHES)
    output:
        f"{data_dir}/done_files/4_concat_maple_output.done"
    resources:
        mem_mb=16000,
        runtime=120 # in minutes
    shell:
        """
        mkdir -p {data_dir}/4
        mkdir -p {out_err_dir}/4
        {pyenv_path} {scripts_dir}/4_concat_maple_output.py \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --n_batch {n_batch} \
            --compress {compress} \
            > {out_err_dir}/4/4_cmo_sm.out 2> {out_err_dir}/4/4_cmo_sm.err
        
        touch {output}
        """

# Step 5a: Process MAPLE placements of the dropout masked dataset
rule process_maple_placements_masked:
    threads: 1
    input:
        f"{data_dir}/done_files/4_concat_maple_output.done"
    output:
        f"{data_dir}/done_files/5_process_maple_placements_masked.done"
    resources:
        mem_mb=32000,
        runtime=240 # in minutes
    shell:
        """
        mkdir -p {data_dir}/5
        mkdir -p {out_err_dir}/5
        {pyenv_path} {scripts_dir}/5_process_maple_placements.py \
            --masked_or_random masked \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --samples_dir {samples_dir} \
            --path_metadata_tsv {metadata} \
            --n_diff_mut {n_diff_mut} \
            --masked_max_dist {masked_max_dist} \
            --masking_ratio {masking_ratio} \
            --het_thr {het_thr} \
            --compress {compress} \
            > {out_err_dir}/5/5_pMSP_m_sm.out 2> {out_err_dir}/5/5_pMSP_m_sm.err
        
        touch {output}
        """

# Step 5b: Process MAPLE placements of the randomly masked dataset
rule process_maple_placements_random:
    threads: 1
    input:
        f"{data_dir}/done_files/4_concat_maple_output.done"
    output:
        f"{data_dir}/done_files/5_process_maple_placements_random.done"
    resources:
        mem_mb=32000,
        runtime=240 # in minutes
    shell:
        """
        mkdir -p {data_dir}/5
        mkdir -p {out_err_dir}/5
        {pyenv_path} {scripts_dir}/5_process_maple_placements.py \
            --masked_or_random random \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --samples_dir {samples_dir} \
            --path_metadata_tsv {metadata} \
            --n_diff_mut {n_diff_mut} \
            --masked_max_dist {masked_max_dist} \
            --masking_ratio {masking_ratio} \
            --het_thr {het_thr} \
            --compress {compress} \
            > {out_err_dir}/5/5_pMSP_r_sm.out 2> {out_err_dir}/5/5_pMSP_r_sm.err
        
        touch {output}
        """

# Step 6a: Find contaminant candidates of the dropout masked dataset
rule find_contaminant_candidates_masked:
    threads: num_cores
    input:
        f"{data_dir}/done_files/5_process_maple_placements_masked.done"
    output:
        f"{data_dir}/done_files/6_find_contaminants_candidates_masked.done"
    resources:
        mem_mb=256000,
        runtime=800 # in minutes
    shell:
        """
        mkdir -p {out_err_dir}/6
        mkdir -p {data_dir}/6
        {pypy_path} {scripts_dir}/6_find_contaminants_candidates.py \
            --masked_or_random masked \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --samples_dir {samples_dir} \
            --path_ref_seq {path_ref_seq} \
            --num_cores {num_cores} \
            --compress {compress} \
            > {out_err_dir}/6/6_fcc_m_sm.out 2> {out_err_dir}/6/6_fcc_m_sm.err
        
        touch {output}
        """

# Step 6b: Find contaminant candidates of the randomly masked dataset
rule find_contaminant_candidates_random:
    threads: num_cores
    input:
        f"{data_dir}/done_files/5_process_maple_placements_random.done"
    output:
        f"{data_dir}/done_files/6_find_contaminants_candidates_random.done"
    resources:
        mem_mb=256000,
        runtime=800 # in minutes
    shell:
        """
        mkdir -p {out_err_dir}/6
        mkdir -p {data_dir}/6
        {pypy_path} {scripts_dir}/6_find_contaminants_candidates.py \
            --masked_or_random random \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --samples_dir {samples_dir} \
            --path_ref_seq {path_ref_seq} \
            --num_cores {num_cores} \
            --compress {compress} \
            > {out_err_dir}/6/6_fcc_r_sm.out 2> {out_err_dir}/6/6_fcc_r_sm.err

        touch {output}
        """

# Step 7a: Apply adapted Eyre et al. model to the dropout masked dataset
rule apply_adapted_eyre_model_masked:
    threads: 1
    input:
        f"{data_dir}/done_files/6_find_contaminants_candidates_masked.done"
    output:
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_masked.done"
    resources:
        mem_mb=128000,
        runtime=800 # in minutes
    shell:
        """
        mkdir -p {out_err_dir}/7
        mkdir -p {data_dir}/7
        {pyenv_path} {scripts_dir}/7_adapted_eyre_model.py \
            --masked_or_random masked \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --compress {compress} \
            > {out_err_dir}/7/7_aem_m_sm.out 2> {out_err_dir}/7/7_aem_m_sm.err
        
        touch {output}
        """

# Step 7b: Apply adapted Eyre et al. model to the randomly masked dataset
rule apply_adapted_eyre_model_random:
    threads: 1
    input:
        f"{data_dir}/done_files/6_find_contaminants_candidates_random.done"
    output:
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_random.done"
    resources:
        mem_mb=128000,
        runtime=800 # in minutes
    shell:
        """
        mkdir -p {out_err_dir}/7
        mkdir -p {data_dir}/7
        {pyenv_path} {scripts_dir}/7_adapted_eyre_model.py \
            --masked_or_random random \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --compress {compress} \
            > {out_err_dir}/7/7_aem_r_sm.out 2> {out_err_dir}/7/7_aem_r_sm.err
        
        touch {output}
        """

# Step 8: Concatenate Eyre model output, plot figures
rule concat_eyre_output:
    threads: num_cores
    input:
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_masked.done",
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_random.done"
    output:
        f"{data_dir}/done_files/8_concat_eyre_output.done"
    resources:
        mem_mb=128000,
        runtime=30 # in minutes
    shell:
        """
        mkdir -p {out_err_dir}/8
        mkdir -p {data_dir}/8
        {pyenv_path} {scripts_dir}/8_concat_eyre_output.py \
            --metadata {metadata} \
            --data_dir {data_dir} \
            --param_term {param_term} \
            --path_ref_seq {path_ref_seq} \
            --compress {compress} \
            > {out_err_dir}/8/8_ceo_sm.out 2> {out_err_dir}/8/8_ceo_sm.err
        
        touch {output}
        """