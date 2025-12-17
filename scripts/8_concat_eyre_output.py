from pathlib import Path
import plotly.express as px
import pandas as pd
from tqdm import tqdm
import argparse
from _aux_functions import compress_file, build_seq_from_ref_and_maple, load_dict_from_tsv, smart_open


parser = argparse.ArgumentParser(
    description="Process the outputs of the Eyre model and print figures."
)

parser.add_argument("--metadata", type=str,
                    help="Path to the file containing the metadata of all the samples of the pipeline.")
parser.add_argument('--data_dir', type=str,
                    help='Directory containing all the files generated during the pipeline.')
parser.add_argument('--param_term', type=str,
                    help='str built from the initial parameters used to recognize the files.')
parser.add_argument('--path_ref_seq', type=str,
                    help='Path to the reference sequence used for the pipeline.')
parser.add_argument('--compress', type=int,
                    help='1 or 0, enables compression.')
args = parser.parse_args()

# Get the arguments
data_dir = args.data_dir
param_term = args.param_term
path_ref_seq = args.path_ref_seq
metadata_tsv = args.metadata
compress = args.compress

def gen_ecdf(df_to_plot, x_col, color_col, x_label, y_label, color_sequence, max_points=500, show_legend=True, upper_clip=None, lower_clip=None):
    """Plot the empirical cumulative distribution function (ECDF) of a given dataframe column.

    Args:
        df_to_plot (pd.DataFrame): The dataframe containing the data to plot.
        x_col (str): The column name to plot on the x-axis.
        color_col (str): The column name to use for coloring the lines.
        title (str): The title of the plot.
        x_label (str): The label for the x-axis.
        y_label (str): The label for the y-axis.
        color_sequence (list): A list of colors to use for the lines.
    Returns:
        fig: The plotly figure object."""
        
    if upper_clip is not None:
        df_to_plot = df_to_plot[df_to_plot[x_col] <= upper_clip]
    if lower_clip is not None:
        df_to_plot = df_to_plot[df_to_plot[x_col] >= lower_clip]

    fig = px.ecdf(
        df_to_plot,
        x=x_col,
        color=color_col,
        labels={x_col: x_label},
        color_discrete_sequence=color_sequence
    )

    fig.update_layout(xaxis_title=x_label, yaxis_title=y_label,
                      width=700, height=600,
                      margin=dict(l=40, r=40, t=10, b=40),
                      font=dict(color="#000000", size=14),
                        showlegend=show_legend
                        )
    
    for i, trace in enumerate(fig.data):
        n = len(trace.x)
        if n > max_points:
            idx = np.linspace(0, n - 1, max_points, dtype=int)
            trace.x = trace.x[idx]
            trace.y = trace.y[idx]
    return fig

def plot_metrics_results(df_concat, data_dir, param_term):
    """Distributions of all metrics
    """
    color_discrete_sequence=['#3B6FB6', '#734595']
    # N_HET_sites
    fig_n_het_sites = gen_ecdf(df_concat, 'n_het_sites', 'type', 'Number of heterozygous sites',
                               'Proportion of samples', color_discrete_sequence, upper_clip=5000, show_legend=True)
    # DIST_DIFF
    fig_dist_diff = gen_ecdf(df_concat, 'dist_diff', 'type', 'Distance difference with unmasked', 'Proportion of samples',
                             color_discrete_sequence, upper_clip=10, show_legend=False)
    # DISTANCE
    fig_distance = gen_ecdf(df_concat, 'distance', 'type', 'Distance to reference genome', 'Proportion of samples',
                            color_discrete_sequence, show_legend=False)
    # SUPPORT
    fig_support = gen_ecdf(df_concat, 'support', 'type', 'Placement support', 'Proportion of samples',
                           color_discrete_sequence, show_legend=False)
    # PROP_GEN_MASKED
    fig_prop_gen_masked = gen_ecdf(df_concat, 'prop_gen_masked', 'type', 'Proportion of genome masked', 'Proportion of samples',
                                 color_discrete_sequence, show_legend=False)
    # PROP_DIST_REDUCED
    fig_prop_dist_reduced = gen_ecdf(df_concat, 'prop_dist_reduced', 'type', 'Proportion of distance reduced', 'Proportion of samples',
                                     color_discrete_sequence, show_legend=False)
    # n_candidates
    fig_n_candidates = gen_ecdf(df_concat, 'n_candidates', 'type', 'Number of placement candidates', 'Proportion of samples',
                                color_discrete_sequence, show_legend=False)
    # score_1
    fig_closeness = gen_ecdf(df_concat, 'score_1', 'type', 'Score 1', 'Proportion of samples',
                             color_discrete_sequence, show_legend=False)
    # score_2
    fig_closeness = gen_ecdf(df_concat, 'score_2', 'type', 'Score 2', 'Proportion of samples',
                             color_discrete_sequence, show_legend=False)
    # closeness_ratio
    fig_closeness_ratio = gen_ecdf(df_concat, 'final_score', 'type', 'Closeness ratio', 'Proportion of samples',
                                color_discrete_sequence, show_legend=False)
    # Assemble all the figs on 3 x 3 subplots
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=3, cols=3,
                        subplot_titles=("A. Number of heterozygous sites",
                                        "B. Distance difference with unmasked",
                                        "C. Number of candidate contaminants",
                                        "D. Distance to the tree after masking",
                                        "E. Support of the new placement",
                                        "F. Closeness to the best candidates",
                                        "G. Proportion of genome masked",
                                        "H. Proportion of distance reduced",
                                        "I. Closeness ratio"),
                        # Reduce titles font size
                        horizontal_spacing=0.1, vertical_spacing=0.1)
    fig.add_trace(fig_n_het_sites.data[0], row=1, col=1)
    fig.add_trace(fig_n_het_sites.data[1], row=1, col=1)
    fig.add_trace(fig_dist_diff.data[0], row=1, col=2)
    fig.add_trace(fig_dist_diff.data[1], row=1, col=2)
    fig.add_trace(fig_n_candidates.data[0], row=1, col=3)
    fig.add_trace(fig_n_candidates.data[1], row=1, col=3)
    
    fig.add_trace(fig_distance.data[0], row=2, col=1)
    fig.add_trace(fig_distance.data[1], row=2, col=1)
    fig.add_trace(fig_support.data[0], row=2, col=2)
    fig.add_trace(fig_support.data[1], row=2, col=2)
    fig.add_trace(fig_closeness.data[0], row=2, col=3)
    fig.add_trace(fig_closeness.data[1], row=2, col=3)
    
    fig.add_trace(fig_prop_gen_masked.data[0], row=3, col=1)
    fig.add_trace(fig_prop_gen_masked.data[1], row=3, col=1)
    fig.add_trace(fig_prop_dist_reduced.data[0], row=3, col=2)
    fig.add_trace(fig_prop_dist_reduced.data[1], row=3, col=2)
    fig.add_trace(fig_closeness_ratio.data[0], row=3, col=3)
    fig.add_trace(fig_closeness_ratio.data[1], row=3, col=3)
    
    fig.update_layout()
    # fig_closeness_ratio should have x_axis between 0 and 1
    fig.update_xaxes(range=[0, 1], row=3, col=3)
    # Cap closeness at 20 (xaxis)
    fig.update_xaxes(range=[0, 20], row=2, col=3)
    # Cap n_candidates at 200k
    fig.update_xaxes(range=[0, 20], row=1, col=3)
    # Remove legend of all but the first subplot
    for i in range(2, 18):
        fig['data'][i]['showlegend'] = False

    # Update the layout
    fig.update_layout(height=1200, width=1000,
                      # Update the legend so that "masked" and "random" appear only once
                      legend=dict(
                        title="Type",
                        x=1.05,
                        y=1,
                        xanchor="left",
                        yanchor="top"
                    ),
                      title_x=0.5,
                      margin=dict(l=40, r=40, t=40, b=40),
                      font=dict(color="#000000", size=12))
    
    fig.update_annotations(font_size=10)
    # Also save as png
    fig.write_image(figs_dir / "ecdf_all_metrics.png", scale=2)

if __name__ == "__main__":
    psmp_path_masked = Path(f"{data_dir}/6/processed_placements_results_with_contaminants_masked_{param_term}.tsv")
    psmp_path_random = Path(f"{data_dir}/6/processed_placements_results_with_contaminants_random_{param_term}.tsv")
    psmp_path_masked_df = pd.read_csv(psmp_path_masked, sep="\t")
    psmp_path_random_df = pd.read_csv(psmp_path_random, sep="\t")
    
    eyre_model_masked_output = Path(f"{data_dir}/7/eyre_model_masked_output.tsv")
    eyre_model_random_output = Path(f"{data_dir}/7/eyre_model_random_output.tsv")
    masked_results_tsv = pd.read_csv(eyre_model_masked_output, sep="\t")
    random_results_tsv = pd.read_csv(eyre_model_random_output, sep="\t")

    # Check the crossing between masked and random samples
    set_cont_samples_masked = set(masked_results_tsv[masked_results_tsv['contaminated'] == 1]['id'].tolist())
    set_cont_samples_random = set(random_results_tsv[random_results_tsv['contaminated'] == 1]['id'].tolist())
    common_contaminated = set_cont_samples_masked.intersection(set_cont_samples_random)
    
    print(f"Number of contaminated samples (masked model): {len(set_cont_samples_masked)} out of {masked_results_tsv.shape[0]}")
    print(f"Number of contaminated samples (random model): {len(set_cont_samples_random)} out of {random_results_tsv.shape[0]}")
    print(f"Number of common contaminated samples between masked and random models: {len(common_contaminated)}")

    # # This commented part generates the FASTA for lineage assessment
    # final_fasta_masked_path = Path(
    #     f"{data_dir}/8/contaminated_samples_masked_{param_term}.fasta"
    # )
    # final_fasta_unmasked_path = Path(
    #     f"{data_dir}/8/contaminated_samples_unmasked_{param_term}.fasta"
    # )
    
    # final_fasta_masked_path.parent.mkdir(parents=True, exist_ok=True)
    # maple_vcf_dict_path = Path(f"{data_dir}/2/dict/maple_vcf_dict_{param_term}_masked.tsv")
    # maple_vcf_dict = load_dict_from_tsv(maple_vcf_dict_path)
    # # load reference sequence
    # with open(path_ref_seq, "r") as f:
    #     ref_seq_lines = f.readlines()
    # ref_seq = "".join([line.strip() for line in ref_seq_lines if not line.startswith(">")]).upper()
    
    # # Output a FASTA file with all the masked and all the random samples that are in the contaminated sets
    # list_masked_seqs = []
    # list_unmasked_seqs = []
    
    # for sample in set_cont_samples_masked:
    #     sample_masked_maple_entry = maple_vcf_dict.get((sample, "masked"), [])
    #     sample_masked_seq = build_seq_from_ref_and_maple(ref_seq, sample_masked_maple_entry)
        
    #     sample_unmasked_maple_entry = maple_vcf_dict.get((sample, "unmasked"), [])
    #     sample_unmasked_seq = build_seq_from_ref_and_maple(ref_seq, sample_unmasked_maple_entry)
        
    #     list_masked_seqs.append((sample, sample_masked_seq))
    #     list_unmasked_seqs.append((sample, sample_unmasked_seq))
    
    # with open(final_fasta_masked_path, "w") as f:
    #     for sample_id, seq in tqdm(list_masked_seqs, desc="Writing masked contaminated FASTA"):
    #         f.write(f">{sample_id}\n")
    #         # Write sequence in lines of 70 characters
    #         for i in range(0, len(seq), 70):
    #             f.write(seq[i:i+70] + "\n")
    
    # with open(final_fasta_unmasked_path, "w") as f:
    #     for sample_id, seq in tqdm(list_unmasked_seqs, desc="Writing unmasked contaminated FASTA"):
    #         f.write(f">{sample_id}\n")
    #         # Write sequence in lines of 70 characters
    #         for i in range(0, len(seq), 70):
    #             f.write(seq[i:i+70] + "\n")
    
    print("Now plotting figures using the metadata of contaminated samples")
    
    # Open metadata tsv and retrieve the country of origin of the contaminated samples to build a figure later
    with smart_open(metadata_tsv, "rt") as f:
        metadata_df = pd.read_csv(f, sep="\t", low_memory=False)
    # Only keep columns Run and Country and rename them
    metadata_df = metadata_df[['Run', 'Country', 'Viridian_pangolin_1.29']]
    contaminated_samples_info = metadata_df[metadata_df['Run'].isin(set_cont_samples_masked)]
    
    # Retrieve counts per country
    country_counts = contaminated_samples_info['Country'].value_counts().reset_index()
    country_counts.columns = ['Country', 'Count']
    # choropleth prevalence map of samples
    fig = px.choropleth(
        country_counts,
        locations="Country",
        locationmode="country names",
        color="Count",
        title="Geographical distribution of contaminated samples (masked model)",
        hover_name="Country",
        hover_data={"Count": True, "Country": True},
        color_continuous_scale=px.colors.sequential.Plasma,
    )
    
    figs_dir = Path(data_dir) / "figs"
    print(f"Generating figures in {figs_dir}")
    fig.write_html(f"{figs_dir}/contaminated_samples_geographical_distribution_masked_{param_term}.html")
    
    # Retrieve counts per lineage
    lineage_counts = contaminated_samples_info['Viridian_pangolin_1.29'].value_counts().reset_index()
    lineage_counts.columns = ['Viridian_pangolin_1.29', 'Count']
    # Remove lineages with counts < 3
    lineage_counts = lineage_counts[lineage_counts['Count'] > 2]
    fig_pie = px.pie(
        lineage_counts,
        names='Viridian_pangolin_1.29',
        values='Count',
        title='Prevalent lineages among the contaminated samples',
        # color='Viridian_pangolin_1.29',
        # hover_name="Viridian_pangolin_1.29",
        # hover_data={"Count": True, "Viridian_pangolin_1.29": True},
    )

    fig_pie.write_html(f"{figs_dir}/contaminated_samples_lineage_distribution_masked_{param_term}.html")


    print(f'The two dataframes have shapes: {masked_results_tsv.shape}, {random_results_tsv.shape}')
    print(f"The Dataframes have columns: {masked_results_tsv.columns.tolist()}")
    print('Concatenating the two dataframes')

    df_concat = pd.concat([masked_results_tsv, random_results_tsv], axis=0)

    # Check weird stuff (negative dist diff and hamming dist for instance)
    print(f"Number of samples with negative dist_diff: {df_concat[df_concat['dist_diff'] < 0].groupby('type').size()}")
    print(f"Number of samples with negative closeness_ratio: {df_concat[df_concat['final_score'] < 0].groupby('type').size()}")
    # Prop dist reduced > 1
    print(f"Number of samples with prop_dist_reduced > 1: {df_concat[df_concat['prop_dist_reduced'] > 1].groupby('type').size()}")

    # FIGS
    plot_metrics_results(df_concat=df_concat,
                         data_dir=data_dir,
                         param_term=param_term)
    
    print(f"Now classifying samples using an Isolation Forest")
    col= ['sample_name', 'unmasked', 'distance', 'n_masked_cons', 'n_masked',
       'n_het_sites', 'prop_gen_masked', 'prop_dist_reduced',
       'unmasked_mutations_to_tree', 'mutations_to_tree',
       'unmasked_placement', 'unmasked_support', 'placement', 'support',
       'mutations_masked', 'type', 'dist_diff', 'masking_ratio']

    col_to_keep = ['sample_name', 'type', 'distance', 'n_het_sites', 'prop_gen_masked',
                   'prop_dist_reduced', 'support', 'dist_diff', 'n_candidates', 'score_1', 'score_2', 'final_score']
    feature_cols = ['distance', 'n_het_sites', 'prop_gen_masked', 'prop_dist_reduced', 'support',
                    'dist_diff', 'n_candidates', 'score_1', 'score_2', 'final_score']
    
    # Feature cols to dtype np.float32
    df_concat[feature_cols] = df_concat[feature_cols].astype(np.float32)
    
    masked_results_tsv[feature_cols] = masked_results_tsv[feature_cols].astype(np.float32)
    random_results_tsv[feature_cols] = random_results_tsv[feature_cols].astype(np.float32)

    clf = IsolationForest(contamination='auto', random_state=9)
    clf.fit(masked_results_tsv[feature_cols])
    masked_results_tsv['scores'] = -clf.decision_function(masked_results_tsv[feature_cols])
    random_results_tsv['scores'] = -clf.decision_function(random_results_tsv[feature_cols])
    
    # The higher scores correspond to the outliers --> samples that are more different from the random masking
    # Plot the scores distribution
    plt.figure(figsize=(8, 4))
    all_scores = np.concatenate([masked_results_tsv['scores'], random_results_tsv['scores']])
    bins = np.linspace(all_scores.min(), all_scores.max(), 30)  # 30 equally spaced bins
    # Sizes of the two datasets
    sns.histplot(masked_results_tsv['scores'], color='#3B6FB6', label='Masked', bins=bins, kde=True)
    sns.histplot(random_results_tsv['scores'], color='#734595', label='Random', bins=bins, kde=True)

    plt.legend()

    plt.xlabel("Anomaly score")
    plt.ylabel("Count")

    plt.savefig(figs_dir / "anomaly_scores_histogram.png", dpi=300)
    
    print("Distribution of anomaly scores saved")

    path_df_masked = f"{data_dir}/8/masked_samples_{param_term}.tsv" # Final OUTPUT
    path_df_random = f"{data_dir}/8/random_samples_{param_term}.tsv"
    masked_results_tsv.to_csv(path_df_masked, sep="\t", index=False)
    random_results_tsv.to_csv(path_df_random, sep="\t", index=False)
    print(f"Masked samples with scores saved to {path_df_masked}")
    print(f"Random samples with scores saved to {path_df_random}")
