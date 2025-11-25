import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scipy.interpolate import make_interp_spline

def pull_from_cluster(target_list, remote_server_name, username, main_server_path, results_folder_name, main_local_path):
    script_name = 'pull_score_files.sh'
    with open(main_local_path+script_name, 'w') as request_file:
        request_file.write("#!/bin/bash\nsftp -q ")
        request_file.write(username+"@"+remote_server_name+" << EOF")

    with open(main_local_path+script_name, 'a') as request_file:
        for target in target_list:
            if '/' in target:
                index_of_target = target.find('/')
                target_folder = target[:index_of_target] + '_' + target[index_of_target + 1:]
            else:
                target_folder = target
            request_file.write('\nget '+main_server_path+target+'/'+results_folder_name+'score*.sc "'+main_local_path+target_folder+'"')
            if not os.path.isdir(main_local_path+target_folder):
                os.mkdir(main_local_path+target_folder)


def compile_score_data_frame(path, docking_flag):
    """
    This function compiles the score files in the specified path into a pandas dataframe.
    Modified from Augustine Arredono

    :param path: (string) the folder that contains the score files (.sc)
    :param docking_flag: (boolean) indicator for whether the score files came from a Relax protocol (False) or
    Docking protocol (True).
    :return pd.concat(score_list): a pandas dataframe of all the score files' info compiled
    """

    filenames = [filename for filename in os.listdir(path) if '.sc' in filename]  # list comprehension to compile .sc
    score_list = []  # an empty variable to compile all the files contents

    for name in filenames:  # Reading all the .sc files
        if docking_flag:
            header_value = 0
        else: # header=1 for relax score files
            header_value = 1
        dataframe = pd.read_csv(path + name, header=header_value, sep=r'\s+')  
        if 'SCORE:' in dataframe.columns:
            del dataframe['SCORE:']
        score_list.append(dataframe)
    final_dataframe = pd.concat(score_list)
    final_dataframe['total_interface_energy'] = np.sum(
        final_dataframe.loc[:, final_dataframe.columns.str.contains('interf_')], axis=1)
    return final_dataframe

def plot_constraint_distribution(dataframe, cst_columns, threshold, subfolder_path, folder_name):
    """
    Plots histogram and cumulative histogram (as subplots) for all constraint columns.
    """
    bin_width = 1
    num_features = len(cst_columns)

    # --- Histogram subplots ---
    fig_hist, axes_hist = plt.subplots(num_features, 1, figsize=(10, 5 * num_features), sharex=False)
    if num_features == 1:
        axes_hist = [axes_hist]
    for ax, feature in zip(axes_hist, cst_columns):
        data = dataframe[feature].dropna().values
        bins = np.arange(min(data), max(data) + bin_width, bin_width)
        ax.hist(data, bins=bins,color='skyblue',label=f'{feature}')
        ax.axvline(threshold, color='black', linestyle='--', label=f'Threshold ({threshold:.2f})')
        ax.set_title(f"Histogram of {feature}")
        ax.set_xlabel(feature)
        ax.set_ylabel("Count")
        ax.set_xlim(left=0)
        ax.legend()
    fig_hist.suptitle(f"{folder_name}: Constraint Histograms")
    fig_hist.tight_layout()
    fig_hist.subplots_adjust(top=0.93)
    fig_hist.savefig(os.path.join(subfolder_path, f"{folder_name}_constraint_histograms.png"), dpi=300)
    plt.close(fig_hist)

    # --- Cumulative histogram subplots ---
    fig_cum, axes_cum = plt.subplots(num_features, 1, figsize=(10, 5 * num_features), sharex=False)
    if num_features == 1:
        axes_cum = [axes_cum]
    for ax, feature in zip(axes_cum, cst_columns):
        data = dataframe[feature].dropna().values
        bins = np.arange(min(data), max(data) + bin_width, bin_width)
        counts, bin_edges = np.histogram(data, bins=bins)
        cumulative = np.cumsum(counts)
        percent_data = cumulative / cumulative[-1] * 100
        x = bin_edges[:-1]
        ax.scatter(x, percent_data, color='skyblue', label=f'{feature} (points)')
        ax.axvline(threshold, color='black', linestyle='--', label=f'Threshold ({threshold:.2f})')
        ax.set_title(f"{feature}")
        ax.set_xlabel(feature)
        ax.set_ylabel("Cumulative Percentage")
        ax.set_ylim(bottom=0)
        ax.legend()
    fig_cum.suptitle(f"{folder_name}: Cumulative Constraint Histograms")
    fig_cum.tight_layout()
    fig_cum.subplots_adjust(top=0.93)
    fig_cum.savefig(os.path.join(subfolder_path, f"{folder_name}_constraint_cumulative_histograms.png"), dpi=300)
    plt.close(fig_cum)

def scatter_hist(x, y, ax, ax_histx, ax_histy, bins=40):
    # Scatter plot
    ax.scatter(x, y, color='skyblue')

    # X histogram
    n_x, bins_x, _ = ax_histx.hist(x, bins=bins, color='skyblue', edgecolor='white')
    peak_idx_x = np.argmax(n_x)
    peak_bin_x = (bins_x[peak_idx_x] + bins_x[peak_idx_x+1]) / 2
    total_upto_peak = n_x[:peak_idx_x + 1].sum()

    # Arrow and annotation pointing on the x-axis of scatter plot (below)
    ax.annotate(
        f"Peak at\n{total_upto_peak/n_x.sum():.2f}",
        xy=(peak_bin_x, ax.get_ylim()[0]),        # Point at x peak, bottom of y-axis
        xytext=(peak_bin_x, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.1),  # offset downward
        arrowprops=dict(facecolor='#985d93', shrink=0.05, width=1, headwidth=6),
        ha='center', va='top', color='#985d93'
    )

    # Y histogram
    n_y, bins_y, _ = ax_histy.hist(y, bins=bins, orientation='horizontal', color='skyblue', edgecolor='white')
    peak_idx_y = np.argmax(n_y)
    peak_bin_y = (bins_y[peak_idx_y] + bins_y[peak_idx_y+1]) / 2

    # Arrow and annotation pointing on the y-axis of scatter plot (left)
    ax.annotate(
        f"Peak at\n{1-peak_bin_y/min(y):.2f}",
        xy=(ax.get_xlim()[0], peak_bin_y),        # Point at bottom of x-axis, y peak
        xytext=(ax.get_xlim()[0] - (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1, peak_bin_y),  # offset leftward
        arrowprops=dict(facecolor='#985d93', shrink=0.05, width=1, headwidth=6),
        ha='right', va='center', color='#985d93'
    )

    # Hide tick labels on marginal axes
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

def generate_diagnostic_plots(full_df, filtering_features, thresholds, subfolder_path, folder_name):
    """
    Generates histograms and cumulative sum plots for given features, uses constraint function for 'cst' in name.
    trying the 2D version
    """

    # Use constraint plotting style for constraint columns (as detected by name)
    if isinstance(filtering_features[0], list) and 'cst' in filtering_features[0][0].lower():
        plot_constraint_distribution(full_df, filtering_features[0], thresholds[0], subfolder_path, folder_name)
    else:
        plot_constraint_distribution(full_df, [filtering_features[0]], thresholds[0], subfolder_path, folder_name)
    
    # 2D scatter + marginal histograms for the other two features
    data = full_df[filtering_features[1:]].dropna().values
    x = data[:, 0]
    y = data[:, 1]

    # Set up grids
    fig = plt.figure(figsize=(8, 8))
    # axes positions in [left, bottom, width, height] in 0-1 figure coords
    ax = fig.add_axes([0.1, 0.1, 0.65, 0.65])
    ax_histx = fig.add_axes([0.1, 0.75, 0.65, 0.2])
    ax_histy = fig.add_axes([0.75, 0.1, 0.2, 0.65])

    scatter_hist(x, y, ax, ax_histx, ax_histy, bins=40)

    ax.set_xlabel(f'First filter: {filtering_features[1]}')
    ax.set_ylabel(f'Second filter: {filtering_features[2]}')
    ax.set_title(f"{folder_name}: Energy filtering quality check")

    # Optionally: add threshold lines
    ax.axvline(thresholds[1], color='black', linestyle='--', label=f'Threshold X ({thresholds[1]:.2f})')
    ax.axhline(thresholds[2], color='black', linestyle='--', label=f'Threshold Y ({thresholds[2]:.2f})')
    ax.legend()

    plt.savefig(os.path.join(subfolder_path, f"{folder_name}_energy_filter_quality_check.png"), dpi=300,bbox_inches='tight' )
    plt.close()

    # Cumulative sum for the other two features
    for feature, threshold in zip(filtering_features[1:], thresholds[1:]):
        data = full_df[feature].dropna().values
        counts, bin_edges = np.histogram(data)
        cumsum_data = np.cumsum(counts)
        percent_data = cumsum_data / cumsum_data[-1] * 100
        x = bin_edges[:-1]
        plt.figure(figsize=(8, 5))
        plt.scatter(x, percent_data, color='skyblue', label='Cumulative Sum')
        plt.axvline(threshold, color='black', linestyle='--', label=f'Threshold ({threshold:.2f})')
        plt.title(f"{folder_name}: Cumulative Percentage of {feature}")
        plt.xlabel(feature)
        plt.ylim(bottom=0)
        plt.ylabel("Cumulative Percentage")
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter())
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(subfolder_path, f"{folder_name}_{feature}_cumsum.png"), dpi=300)
        plt.close()


def filter_rosetta_pdbs(score_dataframe, docking_flag, index_order, all_cst_flag, cst_name, constraint_threshold,
                        interface_energy_threshold_percentage, total_energy_threshold_percentage, iteration_name,
                        display_flag, main_local_path, quality_check_flag):
    """
    Filters Rosetta PDBs based on constraints and energy thresholds.
    """
    # processing parameters
    minimum_pdbs = 0  # Bandaid fix
    column_names = ['total_interface_energy', 'total_score']
    column_threshold_percent = [interface_energy_threshold_percentage, total_energy_threshold_percentage]
    current_dataframe = score_dataframe.copy(deep=True)

    # Reorder columns and thresholds based on index_order
    current_order = ["cst"]+[column_names[index] for index in index_order]
    all_features = [column_names[index] for index in index_order]
    column_threshold_percent = [constraint_threshold]+[column_threshold_percent[index] for index in index_order]

    filter_thresholds = [] # order will be cst and then the other column_names in index_order
    filter_counts = [] # order will be cst and then the other column_names in index_order

    if docking_flag:
        current_dataframe['total_interface_energy'] = current_dataframe.filter(like='interf_').sum(axis=1)
        for index in range(len(current_order)):
            current_feature = current_order[index]
            if index == 0 and all_cst_flag:
                current_threshold = constraint_threshold
                boolean_filter = current_dataframe[cst_name] <= current_threshold
                all_features = [cst_name] + all_features
            elif index == 0 and not all_cst_flag:
                current_threshold = constraint_threshold
                cst_columns = current_dataframe.filter(like=f"_{cst_name}")
                boolean_filter = (cst_columns <= current_threshold).all(axis=1)
                all_features = [cst_columns.columns.tolist()] + all_features
            else:
                threshold_percentage = column_threshold_percent[index]
                current_threshold = current_dataframe[current_feature].nsmallest(
                    round(current_dataframe.shape[0] * threshold_percentage)).max()
 
                boolean_filter = current_dataframe[current_feature] <= current_threshold
            filter_thresholds.append(current_threshold)
            filter_counts.append(boolean_filter.sum())
            
            if filter_counts[index] >= minimum_pdbs:
                current_dataframe = current_dataframe[boolean_filter]
            else:
                print(f"{iteration_name}: Ended filtering early due to low number of PDBs after filtering {current_feature}.")
                break

        if display_flag:
            for index in range(len(current_order)):
                if index == 0 and all_cst_flag:
                    current_feature = 'combined_all_cst'
                elif index == 0 and not all_cst_flag:
                    current_feature = 'individual_cst_columns'
                else:
                    current_feature = current_order[index]
                current_threshold = filter_thresholds[index]
                passed_filter_count = filter_counts[index]
                print(f"{passed_filter_count} passed the {current_feature} filter of {current_threshold}")
        if quality_check_flag:
            generate_diagnostic_plots(score_dataframe,all_features, filter_thresholds, main_local_path, iteration_name)

    return current_dataframe

