import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---- User defined variables: ----
main_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/Training/Jennifer/'
desired_columns = None  # or None for all numeric
group_by_column = 'gene_name'  # categorical grouping column
reference_name = 'WT'  # reference variant name to include in each plot
output_folder = os.path.join(main_path, 'plots')
sort_flag = True
single_file = False  # Set to False when working with multiple CSVs in folder
csv_folder = [os.path.join(main_path,'2025_09_22_D2D_database'),
              os.path.join(main_path,'2025_09_22_D2D_database_complete'),
              os.path.join(main_path,'2025_09_22_D2D_database_complete_previous')]  # Folder(s) containing all CSVs if single_file is False
# ---------------------------------

def extract_base_design(name):
    """Extract base design by removing last character from gene_name"""
    if isinstance(name, str) and len(name) > 1:
        return name[:-1]
    else:
        return name

def plot_box_per_base_design_with_individual_genes(
    df,
    reference_name,
    variant_column='gene_name',
    columns=None,
    output_folder=None,
    prefix=''
):
    if columns is None:
        columns = df.select_dtypes(include='number').columns.tolist()

    # Create a filtered dataframe excluding rows with reference_name in variant_column
    df_no_ref = df[df[variant_column] != reference_name].copy()

    # Extract base design by removing last character on filtered dataframe
    df_no_ref['base_design'] = df_no_ref[variant_column].apply(
        lambda x: x[:-1] if isinstance(x, str) and len(x) > 1 else x
    )

    # Unique base designs excluding the reference's rows
    base_designs = df_no_ref['base_design'].unique()

    # Get rows exactly matching the reference_name (unfiltered)
    ref_df = df[df[variant_column] == reference_name]

    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    for base_design in base_designs:
        plt.figure(figsize=(10, 6))

        # All gene_name rows with this base_design (including reference_name excluded above)
        group_df = df_no_ref[df_no_ref['base_design'] == base_design]
        # Combine these for plotting
        plot_df = pd.concat([group_df, ref_df])
        genes_to_plot = plot_df[variant_column].unique().tolist()
        group_data = []
        group_labels = []
        ref_box_index = None

        for idx, gene in enumerate(genes_to_plot):
            vals = plot_df[plot_df[variant_column] == gene][columns[0]].dropna().values
            group_data.append(vals)
            group_labels.append(gene)
            if gene == reference_name:
                ref_box_index = idx

        if sort_flag:
            medians = [np.median(data) for data in group_data]
            sorted_indices = sorted(range(len(medians)), key=lambda i: medians[i])
            group_data = [group_data[i] for i in sorted_indices]
            group_labels = [group_labels[i] for i in sorted_indices]
            new_group_labels = []
            for label in group_labels:
                if label != reference_name:
                    # Replace with last character
                    new_group_labels.append(label[-1])
                else:
                    # Keep reference name unchanged
                    new_group_labels.append(label)

            if ref_box_index is not None:
                ref_box_index = sorted_indices.index(ref_box_index)

        bp = plt.boxplot(group_data, patch_artist=True, tick_labels=new_group_labels)

        for i, box in enumerate(bp['boxes']):
            color = '#974068' if i == ref_box_index else 'skyblue'
            plt.setp(box, facecolor=color, edgecolor='black')

        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['caps'], color='black')
        plt.setp(bp['medians'], color='black')
        plt.setp(bp['fliers'], markerfacecolor='black', markeredgecolor='black')

        #plt.xticks(rotation=45, ha='right')
        plt.ylim(0)
        plt.xlabel('Gene Name')
        plt.ylabel(columns[0])
        plt.title(f"Base design {base_design}: genes vs {reference_name} - {columns[0]}")
        plt.tight_layout()

        if output_folder:
            filename = f"{prefix}{base_design}_genes_vs_{reference_name}_boxplot.png"
            save_path = os.path.join(output_folder, filename)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

def plot_all_in_one_folder_with_base_designs(
    csv_folder,
    reference_name,
    variant_column='gene_name',
    columns=None,
    output_folder=None
):
    csv_folder = csv_folder[0] if isinstance(csv_folder, list) else csv_folder
    for filename in os.listdir(csv_folder):
        if filename.endswith('.csv'):
            csv_path = os.path.join(csv_folder, filename)
            print(f"Processing {csv_path}...")
            df = pd.read_csv(csv_path)
            prefix = os.path.splitext(filename)[0] + "_"
            plot_box_per_base_design_with_individual_genes(
                df,
                reference_name=reference_name,
                variant_column=variant_column,
                columns=columns,
                output_folder=output_folder,
                prefix=prefix
            )

def plot_combining_folders_with_base_designs(
    csv_folders,          # List of folder paths
    reference_name,
    variant_column='gene_name',
    columns=None,
    output_folder=None
):
    # Step 1: Find all unique CSV filenames across all folders
    unique_filenames = set()
    for folder in csv_folders:
        for f in os.listdir(folder):
            if f.endswith('.csv'):
                unique_filenames.add(f)

    # Step 2-4: For each file, combine and plot
    for filename in unique_filenames:
        csvs_to_combine = []
        for folder in csv_folders:
            path = os.path.join(folder, filename)
            if os.path.isfile(path):
                print(f"Reading {path}")
                csvs_to_combine.append(pd.read_csv(path))

        if csvs_to_combine:
            combined_df = pd.concat(csvs_to_combine, ignore_index=True)

            prefix = os.path.splitext(filename)[0] + "_"
            plot_box_per_base_design_with_individual_genes(
                combined_df,
                reference_name=reference_name,
                variant_column=variant_column,
                columns=columns,
                output_folder=output_folder,
                prefix=prefix
            )
        else:
            print(f"No files to combine for {filename}")


if __name__ == "__main__":
    if single_file:
        df = pd.read_csv(os.path.join(main_path, 'all_scores.csv'))
        plot_box_per_base_design_with_individual_genes(
            df,
            reference_name=reference_name,
            variant_column=group_by_column,
            columns=desired_columns,
            output_folder=output_folder
        )
    elif len(csv_folder) < 2:
        plot_all_in_one_folder_with_base_designs(
            csv_folder=csv_folder,
            reference_name=reference_name,
            variant_column=group_by_column,
            columns=desired_columns,
            output_folder=output_folder
        )
    else:
        plot_combining_folders_with_base_designs(
            csv_folders=csv_folder,          # List of folder paths
            reference_name=reference_name,
            variant_column=group_by_column,
            columns=desired_columns,
            output_folder=output_folder
        )

