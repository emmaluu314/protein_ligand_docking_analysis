import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---- User defined variables: ----
main_path = '/Users/emluu/Documents/Siegel lab/Manuscript work/ALDH/ALDH_phase_2/'
desired_columns = ['duration_hours_single']  # list of columns to plot, or None for all numeric
group_by_column = 'protocol'            # group by protocol
plot_title = "Duration Distributions"
single_file = True
sort_flag = True

# if you want specific labels and file:
xlabel_text = "Protocols"
ylabel_text = "Run time (averaged over replicates) in hours"
file_path = os.path.join(main_path, 'all_timings_master.csv') # if single_file is True
ylim = [0,2]
# ---------------------------------

def plot_box_from_csv(
    csv_path,
    columns=None,
    groupby=None,
    save_path=None,
    output_folder=None,
    title=None,
    sort_by_median=False,
    use_seaborn=True,
    xlabel_title=None,
    ylabel_title=None,
    ylim_setting=None
):
    """
    Plots grouped boxplots from a CSV, optionally using seaborn for clarity.
    """
    df = pd.read_csv(csv_path)

    # Normalize columns input
    if isinstance(columns, str):
        columns = [columns]
    if columns is None:
        columns = df.select_dtypes(include='number').columns.tolist()
    missing_cols = [c for c in columns if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing columns: {missing_cols}. Available: {list(df.columns)}")

    # Determine figure width based on N groups (at least min_width)
    if groupby and groupby in df.columns:
        groups = df[groupby].unique()
        num_elements = len(groups)
    else:
        num_elements = len(columns)
    min_width = 5
    width = max(min_width, min_width * max(1, (num_elements // 10)))
    height = 6

    plt.figure(figsize=(width, height))

    col = columns[0]   # Only support one column for y for grouped plot

    if groupby and groupby in df.columns:
        if use_seaborn:
            try:
                import seaborn as sns
            except ImportError:
                raise ImportError("Seaborn is not installed. Set use_seaborn=False to use only Matplotlib.")
            data = df[[groupby, col]].dropna()
            order = None
            if sort_by_median:
                med = data.groupby(groupby)[col].median().sort_values()
                order = med.index.tolist()
            sns.boxplot(x=groupby, y=col, data=data, order=order, showfliers=True)
            if len(med) > 5:
                plt.xticks(rotation=45, ha='right')
        else:
            group_labels = []
            group_data = []
            for g in df[groupby].unique():
                vals = df[df[groupby] == g][col].dropna().values
                group_labels.append(str(g))
                group_data.append(vals)
            if sort_by_median:
                medians = [np.median(data) for data in group_data]
                sorted_idx = np.argsort(medians)
                group_labels = [group_labels[i] for i in sorted_idx]
                group_data = [group_data[i] for i in sorted_idx]
            bp = plt.boxplot(group_data, labels=group_labels, patch_artist=True)
            plt.setp(bp['boxes'], facecolor='skyblue', edgecolor='black')
            plt.setp(bp['whiskers'], color='black')
            plt.setp(bp['caps'], color='black')
            plt.setp(bp['medians'], color='black')
            plt.setp(bp['fliers'], markerfacecolor='black', markeredgecolor='black')
            if len(medians) > 5:
                plt.xticks(rotation=45, ha='right')
        plt.xlabel(xlabel_title if xlabel_title else f"{groupby}")

    else:
        df[columns].boxplot()
        plt.xlabel(xlabel_title if xlabel_title else "Columns")

    plt.title(title if title else f"Boxplot of {os.path.basename(csv_path)}")
    plt.ylabel(ylabel_title if ylabel_title else col)
    if ylim_setting:
        plt.ylim(ylim_setting)
    plt.tight_layout()

    # Output logic
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
        filename = os.path.splitext(os.path.basename(csv_path))[0] + "_boxplot.png"
        save_path = os.path.join(output_folder, filename)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot to {save_path}")
    else:
        plt.show()

def plot_multiple_csv_boxplots(
    csv_folder,
    columns=None,
    groupby=None,
    output_folder=None,
    save_folder=None,
    flag=False,
    use_seaborn=True,
    xlabel_title=None,
    ylabel_title=None,
    ylim_setting=None
):
    """
    Plots boxplots for all CSV files in a folder, optionally grouped and using seaborn.
    """
    save_dir = save_folder or output_folder

    for filename in os.listdir(csv_folder):
        if filename.endswith('.csv'):
            csv_path = os.path.join(csv_folder, filename)
            save_path = None
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
                plot_name = f"{os.path.splitext(filename)[0]}_boxplot.png"
                save_path = os.path.join(save_dir, plot_name)
            plot_box_from_csv(
                csv_path,
                columns=columns,
                groupby=groupby,
                save_path=save_path,
                title=f"Boxplot of {filename}",
                sort_by_median=flag,
                use_seaborn=use_seaborn, 
                xlabel_title=xlabel_title,
                ylabel_title=ylabel_title
            )

# --------------------- Script Entry ------------------------
if __name__ == "__main__":
    if single_file:
        plot_box_from_csv(
            file_path,
            columns=desired_columns,
            groupby=group_by_column,
            output_folder=main_path,
            title=plot_title,
            use_seaborn=True,  # or False if you want only Matplotlib
            sort_by_median=sort_flag,
            xlabel_title=xlabel_text,
            ylabel_title=ylabel_text,
            ylim_setting=ylim
        )
    else:
        plot_multiple_csv_boxplots(
            main_path,
            columns=desired_columns,
            groupby=group_by_column,
            save_folder=os.path.join(main_path, 'plots'),
            flag=sort_flag,
            use_seaborn=True,  # or False for Matplotlib,
            xlabel_title=xlabel_text,
            ylabel_title=ylabel_text,
            ylim_setting=ylim
        )
