from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import argparse
from _aux_functions import compress_file, smart_open
from sklearn.ensemble import IsolationForest

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


def plot_metrics_results(df_concat, data_dir, param_term):
    """Mix of Violin plots (continuous) and Histograms (discrete/binary)"""
    
    color_discrete_sequence = ['#3B6FB6', '#734595']
    
    fig = make_subplots(rows=2, cols=3,
                        subplot_titles=("Number of heterozygous sites",
                                        "Masking ratio",
                                        "Placement support",
                                        "Unique candidate contaminants",
                                        "Closeness score",
                                        "Posterior probability of contamination"),
                        horizontal_spacing=0.08, vertical_spacing=0.15)

    def add_plot(plot_type, target_col, row, col, upper_clip=None, lower_clip=None, nbinsx=None):
        df_clean = df_concat.copy()
        df_clean[target_col] = pd.to_numeric(df_clean[target_col], errors='coerce')
        
        if upper_clip is not None:
            df_clean = df_clean[df_clean[target_col] <= upper_clip]
        if lower_clip is not None:
            df_clean = df_clean[df_clean[target_col] >= lower_clip]

        types = df_clean['type'].dropna().unique()
        
        for i, t in enumerate(sorted(types, reverse=False)):
            subset = df_clean[df_clean['type'] == t][target_col].dropna()
            
            show_legend = True if (row == 1 and col == 1) else False
            color = color_discrete_sequence[i % len(color_discrete_sequence)]
            
            if plot_type == 'violin':
                fig.add_trace(
                    go.Violin(
                        x=[t.capitalize()] * len(subset), # Category on X
                        y=subset,                         # Metric on Y
                        name=t.capitalize(),
                        marker_color=color,
                        box_visible=True,
                        meanline_visible=True,
                        points=False,
                        opacity=0.8,
                        line_width=1.5,
                        width=0.8,         # CRITICAL FIX: Forces the violins to spread out horizontally
                        spanmode='hard',   # Prevents phantom tails extending past your actual data
                        scalemode='width', # Forces both Masked and Random to take up equal visual width
                        legendgroup=t,
                        showlegend=show_legend
                    ),
                    row=row, col=col
                )
            elif plot_type == 'histogram':
                fig.add_trace(
                    go.Histogram(
                        x=subset,                         # Metric on X
                        name=t.capitalize(),
                        marker_color=color,
                        histnorm='probability',           # Probability on Y
                        opacity=0.8,
                        nbinsx=nbinsx,
                        legendgroup=t,
                        showlegend=show_legend
                    ),
                    row=row, col=col
                )

    # --- ROW 1 ---
    # Because n_het_sites is heavily skewed integers, a Histogram is much better!
    add_plot('violin', 'n_het_sites', 1, 1, upper_clip=50, nbinsx=25)
    
    # These are continuous, so Violins are perfect
    add_plot('violin', 'masking_ratio', 1, 2, upper_clip=50)
    add_plot('violin', 'support', 1, 3)
    
    # --- ROW 2 ---
    add_plot('histogram', 'n_unique_candidates', 2, 1, upper_clip=5, lower_clip=0, nbinsx=6)
    add_plot('violin', 'final_score', 2, 2)
    # add_plot('violin', 'n_mutations_removed_by_masking', 2, 2, upper_clip=25) # Difference between the number of mutations to the tree before and after masking
    
    # Discrete/Binary metrics -> Histograms!
    # add_plot('violin', 'n_mutations_removed_by_masking', 2, 3, lower_clip=0, upper_clip=25) # Difference between the number of mutations to the tree before and after masking
    add_plot('histogram', 'contaminated', 2, 3, nbinsx=2)

    # --- AXES FORMATTING ---
    # For Violins, format the Y-AXIS (Metric is on Y)
    # fig.update_yaxes(range=[0, 20], row=1, col=2)  # masking_ratio
    fig.update_yaxes(range=[0, 1], row=1, col=3)   # support
    fig.update_yaxes(range=[0, 1], row=2, col=2)   # final_score
    # fig.update_yaxes(range=[0, 25], row=2, col=2) # n_mutations_removed_by_masking
    
    # For Histograms, format the X-AXIS (Metric is on X)
    fig.update_xaxes(range=[-0.5, 5.5], tickvals=[0, 1, 2, 3, 4, 5], row=2, col=1) # n_unique_candidates
    fig.update_xaxes(range=[-0.5, 1.5], tickvals=[0, 1], row=2, col=3)             # contaminated
    
    # Add "Probability" label to the Y-axis of the Histograms
    fig.update_yaxes(title_text="Probability", row=1, col=1)
    fig.update_yaxes(title_text="Probability", row=2, col=2)
    fig.update_yaxes(title_text="Probability", row=2, col=3)

    # Global layout settings
    fig.update_layout(
        height=800, width=1200, 
        violinmode='group',
        barmode='group',
        bargap=0.1,
        # title_text="Metric Distributions (Masked vs Random)",
        title_x=0.5,
        legend=dict(title="Type", x=1.02, y=1, xanchor="left", yanchor="top"),
        margin=dict(l=40, r=40, t=60, b=40),
        font=dict(color="#000000", size=12),
        template="plotly_white"
    )

    fig.update_annotations(
        font_size=14,
        # Add vertical spacing to subplot titles to prevent overlap with the violins/histograms
        yshift=10,
    )

    # Save the figure
    figs_dir = Path(data_dir) / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)
    
    fig.write_image(figs_dir / "mixed_selected_metrics.png", scale=2)
    fig.write_html(figs_dir / "mixed_selected_metrics.html")

if __name__ == "__main__":
    # CONCATENATE RESULTS FROM SCRIPT 6 AND 7, PLOT FIGURES AND OUTPUT FINAL RESULTS
    psmp_path_masked = Path(f"{data_dir}/6/processed_placements_results_with_contaminants_masked_{param_term}.tsv")
    psmp_path_random = Path(f"{data_dir}/6/processed_placements_results_with_contaminants_random_{param_term}.tsv")
    psmp_masked_df = pd.read_csv(psmp_path_masked, sep="\t")
    psmp_random_df = pd.read_csv(psmp_path_random, sep="\t")
    
    eyre_model_masked_output = Path(f"{data_dir}/7/eyre_model_masked_output.tsv")
    eyre_model_random_output = Path(f"{data_dir}/7/eyre_model_random_output.tsv")
    masked_results_tsv = pd.read_csv(eyre_model_masked_output, sep="\t")
    random_results_tsv = pd.read_csv(eyre_model_random_output, sep="\t")
    # Some samples in the psmp results are missing in the eyre results because we couldn't find them a contaminant to test.
    # Add those in the eyre df, filling both columns with "No contaminant found"
    extra_cols = ['haplotype_pair:p_pair:ML_mu:ML_epsilon', 'contaminated']

    merged_masked_df = psmp_masked_df.merge(
        masked_results_tsv,
        on="sample_name",
        how="left"
    )
    # Fill missing entries for the extra columns
    merged_masked_df['haplotype_pair:p_pair:ML_mu:ML_epsilon'] = merged_masked_df['haplotype_pair:p_pair:ML_mu:ML_epsilon'].fillna("No contaminant found")
    merged_masked_df['contaminated'] = merged_masked_df['contaminated'].fillna(0)

    merged_random_df = psmp_random_df.merge(
        random_results_tsv,
        on="sample_name",
        how="left"
    )
    # Fill missing entries for the extra columns
    merged_random_df[extra_cols] = merged_random_df[extra_cols].fillna("No contaminant found")

    # Check the crossing between masked and random samples
    set_cont_samples_masked = set(merged_masked_df[merged_masked_df['contaminated'] == 1]['sample_name'].tolist())
    set_cont_samples_random = set(merged_random_df[merged_random_df['contaminated'] == 1]['sample_name'].tolist())
    common_contaminated = set_cont_samples_masked.intersection(set_cont_samples_random)
    
    print(f"Number of contaminated samples (masked model): {len(set_cont_samples_masked)} out of {merged_masked_df.shape[0]}")
    print(f"Number of contaminated samples (random model): {len(set_cont_samples_random)} out of {merged_random_df.shape[0]}")
    print(f"Number of common contaminated samples between masked and random models: {len(common_contaminated)}")
    
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
    fig.write_image(f"{figs_dir}/contaminated_samples_geographical_distribution_masked_{param_term}.png")
    
    # Retrieve counts per lineage
    lineage_counts = contaminated_samples_info['Viridian_pangolin_1.29'].value_counts().reset_index()
    lineage_counts.columns = ['Viridian_pangolin_1.29', 'Count']
    # Remove lineages with counts < 5
    lineage_counts = lineage_counts[lineage_counts['Count'] > 4]
    fig_pie = px.pie(
        lineage_counts,
        names='Viridian_pangolin_1.29',
        values='Count',
        title='Prevalent lineages among the contaminated samples',
        # color='Viridian_pangolin_1.29',
        # hover_name="Viridian_pangolin_1.29",
        # hover_data={"Count": True, "Viridian_pangolin_1.29": True},
    )

    fig_pie.write_image(f"{figs_dir}/contaminated_samples_lineage_distribution_masked_{param_term}.png")

    print(f'The two dataframes have shapes: {merged_masked_df.shape}, {merged_random_df.shape}')
    print(f"The Dataframes have columns: {merged_masked_df.columns.tolist()}")
    print('Concatenating the two dataframes')

    df_concat = pd.concat([merged_masked_df, merged_random_df], axis=0)

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
    
    merged_masked_df[feature_cols] = merged_masked_df[feature_cols].astype(np.float32)
    merged_random_df[feature_cols] = merged_random_df[feature_cols].astype(np.float32)

    clf = IsolationForest(contamination='auto', random_state=9)
    clf.fit(merged_masked_df[feature_cols])
    merged_masked_df['scores'] = -clf.decision_function(merged_masked_df[feature_cols])
    merged_random_df['scores'] = -clf.decision_function(merged_random_df[feature_cols])
    
    # The higher scores correspond to the outliers --> samples that are more different from the random masking
    # Plot the scores distribution
    plt.figure(figsize=(8, 4))
    all_scores = np.concatenate([merged_masked_df['scores'], merged_random_df['scores']])
    bins = np.linspace(all_scores.min(), all_scores.max(), 30)  # 30 equally spaced bins
    # Sizes of the two datasets
    sns.histplot(merged_masked_df['scores'], color='#3B6FB6', label='Masked', bins=bins, kde=True)
    sns.histplot(merged_random_df['scores'], color='#734595', label='Random', bins=bins, kde=True)

    plt.legend()

    plt.xlabel("Anomaly score")
    plt.ylabel("Count")

    plt.savefig(figs_dir / "anomaly_scores_histogram.png", dpi=300)
    
    print("Distribution of anomaly scores saved")

    path_df_masked = f"{data_dir}/8/masked_samples_{param_term}.tsv"
    path_df_random = f"{data_dir}/8/random_samples_{param_term}.tsv"
    merged_masked_df.to_csv(path_df_masked, sep="\t", index=False)
    merged_random_df.to_csv(path_df_random, sep="\t", index=False)
    print(f"Masked samples with scores saved to {path_df_masked}")
    print(f"Random samples with scores saved to {path_df_random}")
