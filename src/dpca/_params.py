
# General
n_batch = 512
pyenv_path = "/homes/anoufa/.pyenv/versions/3.11.6/bin/python3"
pypy_path = "/nfs/research/goldman/anoufa/.venv/bin/pypy3.10"
MAPLE_path = "/nfs/research/goldman/anoufa/src/dpca/MAPLEv0.7.5.py"
param_path = '/nfs/research/goldman/anoufa/src/dpca/_params.py'
samples_dir = '/nfs/research/goldman/anoufa/pipeline/run2viridian_dir.tsv.xz'
num_cores = 8

# 1_gen_maple_file parameters
masking_method="thr"
# thr or hmm 
het_thr = 0.1
# used to filter what's considered a het site
prop_under_depth_thr = 1
# unused
typical_depth_thr = 0
# unused
n_masked_thr = 30000
# unused

depth_thr = 0.1
# MASKING THR, works ~the same with 0.15 and 0.2
max_n_het_sites = 3
# MAX NUMBER OF HET SITES ALLOWED TO CONSIDER A SAMPLE CLEAN ENOUGH TO BE PART OF THE CLEAN TREE (and excluded from the analysis)


param_term = f"{depth_thr}_{max_n_het_sites}"
path_ref_seq = '/nfs/research/goldman/anoufa/data/MAPLE_input/maple_ref_lower.fasta'

# 3_maple_sample_placement parameters
inputTree=f"/nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_aug{param_term}_tree.tree"
inputRates=f"/nfs/research/goldman/anoufa/data/MAPLE_output/clean_tree/clean_aug{param_term}_subs.txt"

# 5_handle_sample_placements parameters
n_diff_mut = 2
masked_max_dist = 5
masking_ratio = 0
min_support = 0
processed_placement_masked = '/nfs/research/goldman/anoufa/pipeline/5/processed_placements_results_masked_0.1_3.1.tsv'
processed_placement_random = '/nfs/research/goldman/anoufa/pipeline/5/processed_placements_results_random_0.1_3.1.tsv'
processed_placement_random_with_contaminants = '/nfs/research/goldman/anoufa/pipeline/5/processed_placements_results_with_contaminants_random_0.1_3.1.tsv'
processed_placement_masked_with_contaminants = '/nfs/research/goldman/anoufa/pipeline/5/processed_placements_results_with_contaminants_masked_0.1_3.1.tsv'
masked_with_scores = '/nfs/research/goldman/anoufa/pipeline/7/masked_samples_0.1_3.tsv'
random_with_scores = '/nfs/research/goldman/anoufa/pipeline/7/random_samples_0.1_3.tsv'
_startprob = [1.0, 0.0, 0.0]
_transmat = [[0.9854401700926695, 0.012928746958638939, 0.001631082948691542], [0.020210870499876913, 0.979535636171372, 0.00025349332875110813], [0.005949357686065885, 0.0005483680201259543, 0.9935022742938081]]
_means = [1205.0, 951.0, 106.0]
_covars = [0.000625, 0.000625, 0.000625]


is_cont_dict_masked = '/nfs/research/goldman/anoufa/data/MAPLE_output/final_filter/is_cont_dict_masked.pickle'
is_cont_dict_random = '/nfs/research/goldman/anoufa/data/MAPLE_output/final_filter/is_cont_dict_random.pickle'
het_sites_dict_masked = '/nfs/research/goldman/anoufa/pipeline/5/het_sites_dict_masked_0.1_3.pickle'
het_sites_dict_random = '/nfs/research/goldman/anoufa/pipeline/5/het_sites_dict_random_0.1_3.pickle'
data_dir = '/nfs/research/goldman/anoufa/pipeline'
maple_alignment_unmasked_file_path = '/nfs/research/goldman/anoufa/pipeline/2/alignment_files/unmasked_alignment_0.1_3.maple'
maple_alignment_random_file_path = '/nfs/research/goldman/anoufa/pipeline/2/alignment_files/random_alignment_0.1_3.maple'
maple_alignment_masked_file_path = '/nfs/research/goldman/anoufa/pipeline/2/alignment_files/masked_alignment_0.1_3.maple'
final_clean_tree_path = '/nfs/research/goldman/anoufa/pipeline/2/alignment_files/clean_tree_alignment_file_0.1_3.maple'
