from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import _params
import numpy as np
from sklearn.ensemble import IsolationForest
from _aux_functions import update_params_file, generate_sh_param_file
from tqdm import tqdm


n_batch = _params.n_batch
het_thr = _params.het_thr
depth_thr = _params.depth_thr
prop_under_depth_thr = _params.prop_under_depth_thr
typical_depth_thr = _params.typical_depth_thr
n_masked_thr = _params.n_masked_thr
path_ref_seq = _params.path_ref_seq
inputTree = _params.inputTree
inputRates = _params.inputRates
n_diff_mut = _params.n_diff_mut
masked_max_dist = _params.masked_max_dist
masking_ratio = _params.masking_ratio
psmp_path_random = _params.processed_placement_random_with_contaminants
psmp_path_masked = _params.processed_placement_masked_with_contaminants
param_term = _params.param_term
param_path = _params.param_path
data_dir = _params.data_dir


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


if __name__ == "__main__":


    df_psmp = pd.read_csv(psmp_path_masked, sep="\t")
    df_psmp_random = pd.read_csv(psmp_path_random, sep="\t")
    
    print(f'The two dataframes have shapes: {df_psmp.shape}, {df_psmp_random.shape}')
    
    print(f"The Dataframes have columns: {df_psmp.columns.tolist()}")
    print('Concatenating the two dataframes')

    df_concat = pd.concat([df_psmp, df_psmp_random], axis=0)

    # Check weird stuff (negative dist diff and hamming dist for instance)
    print(f"Number of samples with negative dist_diff: {df_concat[df_concat['dist_diff'] < 0].groupby('type').size()}")
    print(f"Number of samples with negative closeness_ratio: {df_concat[df_concat['closeness_ratio'] < 0].groupby('type').size()}")
    # Prop dist reduced > 1
    print(f"Number of samples with prop_dist_reduced > 1: {df_concat[df_concat['prop_dist_reduced'] > 1].groupby('type').size()}")

    # FIGS

    figs_folder = Path(f"{data_dir}/7/{param_term}")
    figs_folder.mkdir(parents=True, exist_ok=True)
    print(f"Generating figures in {figs_folder}")
    
    # ECDF of N_HET_SITES, DIST_DIFF, DISTANCE, SUPPORT, MASKING_RATIO, PROP_DIST_REDUCED
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
    # closeness
    fig_closeness = gen_ecdf(df_concat, 'closeness', 'type', 'Closeness', 'Proportion of samples',
                             color_discrete_sequence, show_legend=False)
    # closeness_ratio
    fig_closeness_ratio = gen_ecdf(df_concat, 'closeness_ratio', 'type', 'Closeness ratio', 'Proportion of samples',
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
    fig.write_html(figs_folder / "ecdf_all_metrics.html")
    # Also save as png
    fig.write_image(figs_folder / "ecdf_all_metrics.png", scale=2)

    print(f"Figures saved in {figs_folder}")
    
    print(f"Now classifying samples using an Isolation Forest")
    
    col= ['sample_name', 'unmasked', 'distance', 'n_masked_cons', 'n_masked',
       'n_het_sites', 'prop_gen_masked', 'prop_dist_reduced',
       'unmasked_mutations_to_tree', 'mutations_to_tree',
       'unmasked_placement', 'unmasked_support', 'placement', 'support',
       'mutations_masked', 'type', 'dist_diff', 'masking_ratio']

    col_to_keep = ['sample_name', 'type', 'distance', 'n_het_sites', 'prop_gen_masked',
                   'prop_dist_reduced', 'support', 'dist_diff', 'n_candidates', 'closeness', 'closeness_ratio']
    feature_cols = ['distance', 'n_het_sites', 'prop_gen_masked', 'prop_dist_reduced', 'support',
                    'dist_diff', 'n_candidates', 'closeness', 'closeness_ratio']
    
    # Feature cols to dtype np.float32
    df_concat[feature_cols] = df_concat[feature_cols].astype(np.float32)
    
    df_psmp[feature_cols] = df_psmp[feature_cols].astype(np.float32)
    df_psmp_random[feature_cols] = df_psmp_random[feature_cols].astype(np.float32)

    clf = IsolationForest(contamination='auto', random_state=7)
    clf.fit(df_psmp[feature_cols])
    df_psmp['scores'] = -clf.decision_function(df_psmp[feature_cols])
    df_psmp_random['scores'] = -clf.decision_function(df_psmp_random[feature_cols])
    
    # The higher scores correspond to the outliers --> samples that are more different from the random masking
    # Plot the scores distribution
    plt.figure(figsize=(8, 4))
    all_scores = np.concatenate([df_psmp['scores'], df_psmp_random['scores']])
    bins = np.linspace(all_scores.min(), all_scores.max(), 30)  # 30 equally spaced bins
    # Sizes of the two datasets
    sns.histplot(df_psmp['scores'], color='#3B6FB6', label='Masked', bins=bins, kde=True)
    sns.histplot(df_psmp_random['scores'], color='#734595', label='Random', bins=bins, kde=True)

    plt.legend()

    plt.xlabel("Anomaly score")
    plt.ylabel("Count")

    plt.savefig(figs_folder / "anomaly_scores_histogram.png", dpi=300)
    
    print("Distribution of anomaly scores saved")

    path_df_masked = f"{data_dir}/7/masked_samples_{param_term}.tsv"
    path_df_random = f"{data_dir}/7/random_samples_{param_term}.tsv"
    df_psmp.to_csv(path_df_masked, sep="\t", index=False)
    df_psmp_random.to_csv(path_df_random, sep="\t", index=False)
    print(f"Masked samples with scores saved to {path_df_masked}")
    print(f"Random samples with scores saved to {path_df_random}")
    
    updates = {
        'masked_with_scores': path_df_masked,
        'random_with_scores': path_df_random,
    }

    update_params_file(f"{param_path}", updates)
    generate_sh_param_file()
    print("Updated _params.py with new file paths")

    done_file_path = Path(f"{data_dir}/done_files/7_figs.done")
    done_file_path.parent.mkdir(parents=True, exist_ok=True)
    done_file_path.touch()
    