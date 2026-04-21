########################################
# PhyCD Snakemake pipeline — EXAMPLE VERSION WITH 100 SAMPLES
########################################

# --- Adjustable parameters (via --config or defaults below) ---
n_batch = int(config.get("n_batch", 1)) # Number of jobs to run in parallel
BATCHES = list(range(n_batch))
num_cores = 1 # Number of cores used when multiprocessing within a job

depth_thr         = float(config.get("depth_thr", 0.05) or 0.05)
max_dropout_masked = int(config.get("max_dropout_masked", 0) or 0)
het_thr           = float(config.get("het_thr", 0.1) or 0.1)
max_n_het_sites   = int(config.get("max_n_het_sites", 3) or 3)

param_term = f"{int(depth_thr*100)}_{max_dropout_masked}_{int(het_thr*100)}_{max_n_het_sites}"

print("Parameters received:", flush=True)
print(f"  n_batch:            {n_batch}",            flush=True)
print(f"  depth_thr:          {depth_thr}",          flush=True)
print(f"  max_dropout_masked: {max_dropout_masked}",  flush=True)
print(f"  het_thr:            {het_thr}",             flush=True)
print(f"  max_n_het_sites:    {max_n_het_sites}",     flush=True)

# --- 5_process_sample_placements default parameters ---
n_diff_mut    = 2
masked_max_dist = 5
masking_ratio = 0

# --- Paths: edit these for your local setup ---
root_dir    = "/export/home1/users/mpath/oanoufa/PhyCD/PhyCD"
data_dir    = f"{root_dir}/pipeline/{param_term}/data"
scripts_dir = f"{root_dir}/scripts"
out_err_dir = f"{root_dir}/pipeline/{param_term}/out_err"
samples_dir = f"{root_dir}/input_data/run2viridian_dir_example.tsv"
metadata    = f"{root_dir}/input_data/viridian_example_samples_metadata.tsv"
path_ref_seq = f"{root_dir}/input_data/maple_ref_lower.fasta"

# Interpreter/tool paths
# Having a separate pypy env makes the pipeline run faster but is not necessary
pyenv_path = f"{root_dir}/.phycd_venv/bin/python"
pypy_path  = f"{root_dir}/.phycd_venv/bin/python"
maple_path = f"{scripts_dir}/MAPLEv0.7.5.py"

compress = 1


# --- Rules (identical logic to HPC version, resources reduced) ---

rule all:
    input:
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


rule gen_maple_file:
    threads: 1
    output:
        f"{data_dir}/done_files/1/1_maple_alignment_batch_{{batch}}.done"
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


rule process_gmf_output:
    threads: 1  # no multiprocessing locally
    input:
        expand(f"{data_dir}/done_files/1/1_maple_alignment_batch_{{batch}}.done", batch=BATCHES)
    output:
        f"{data_dir}/done_files/2_process_gmf_output.done"
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


rule generate_clean_tree:
    threads: 1
    input:
        f"{data_dir}/done_files/2_process_gmf_output.done"
    output:
        f"{data_dir}/done_files/2_generate_clean_tree.done"
    shell:
        """
        if [ ! -f "{data_dir}/clean_tree/clean_{param_term}_tree.tree" ]; then
            mkdir -p {data_dir}/clean_tree
            {pypy_path} {maple_path} \
                --rateVariation \
                --model UNREST \
                --saveInitialTreeEvery 200000 \
                --input {data_dir}/2/alignment_files/clean_tree_alignment_file_{param_term}.maple \
                --output {data_dir}/clean_tree/clean_{param_term} \
                --numCores 1 \
                --overwrite \
                --numTopologyImprovements 1 \
                > {out_err_dir}/2/2_gct_sm.out 2> {out_err_dir}/2/2_gct_sm.err
            touch {output}
        else
            echo "Clean tree already exists. Skipping."
            touch {output}
        fi
        """


rule generate_taxonium_tree:
    threads: 1
    input:
        f"{data_dir}/done_files/2_generate_clean_tree.done"
    output:
        f"{data_dir}/done_files/2_newick_to_taxonium.done"
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
            echo "JSONL tree already exists. Skipping."
            touch {output}
        fi
        """


rule maple_sample_placement:
    threads: 1
    input:
        f"{data_dir}/done_files/2_newick_to_taxonium.done"
    output:
        f"{data_dir}/done_files/3/3_maple_sample_placement_{{batch}}.done"
    shell:
        """
        mkdir -p {data_dir}/3
        mkdir -p {out_err_dir}/3
        mkdir -p {data_dir}/done_files/3
        inputTree="{data_dir}/clean_tree/clean_{param_term}_tree.tree"
        inputRates="{data_dir}/clean_tree/clean_{param_term}_subs.txt"
        inputFile="{data_dir}/1/maple_alignment_batch{wildcards.batch}_{param_term}.maple"
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


rule concat_maple_output:
    threads: 1
    input:
        expand(f"{data_dir}/done_files/3/3_maple_sample_placement_{{batch}}.done", batch=BATCHES)
    output:
        f"{data_dir}/done_files/4_concat_maple_output.done"
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


rule process_maple_placements_masked:
    threads: 1
    input:
        f"{data_dir}/done_files/4_concat_maple_output.done"
    output:
        f"{data_dir}/done_files/5_process_maple_placements_masked.done"
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


rule process_maple_placements_random:
    threads: 1
    input:
        f"{data_dir}/done_files/4_concat_maple_output.done"
    output:
        f"{data_dir}/done_files/5_process_maple_placements_random.done"
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


rule find_contaminant_candidates_masked:
    threads: 1
    input:
        f"{data_dir}/done_files/5_process_maple_placements_masked.done"
    output:
        f"{data_dir}/done_files/6_find_contaminants_candidates_masked.done"
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


rule find_contaminant_candidates_random:
    threads: 1
    input:
        f"{data_dir}/done_files/5_process_maple_placements_random.done"
    output:
        f"{data_dir}/done_files/6_find_contaminants_candidates_random.done"
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


rule apply_adapted_eyre_model_masked:
    threads: 1
    input:
        f"{data_dir}/done_files/6_find_contaminants_candidates_masked.done"
    output:
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_masked.done"
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


rule apply_adapted_eyre_model_random:
    threads: 1
    input:
        f"{data_dir}/done_files/6_find_contaminants_candidates_random.done"
    output:
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_random.done"
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


rule concat_eyre_output:
    threads: 1
    input:
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_masked.done",
        f"{data_dir}/done_files/7_apply_adapted_eyre_model_random.done"
    output:
        f"{data_dir}/done_files/8_concat_eyre_output.done"
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