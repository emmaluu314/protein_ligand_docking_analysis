import os
import re
import math
from collections import defaultdict, Counter
import pandas as pd
import matplotlib.pyplot as plt

# ---- User defined variables ----
main_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/LDH/'
csv_folders = [
    os.path.join(main_path, '2025_09_01_16NAD_mutants'),
    os.path.join(main_path, '2025_09_08_16NAD_mutants'),
    os.path.join(main_path, '2025_09_18_16NAD_mutants')
]
output_folder = os.path.join(main_path, 'plots_residue_frequency')
# -------------------------------

def get_residue_chain_presence(contact_list):
    """Return a dict {(resi, chain): 1} for each contact in the list."""
    if not isinstance(contact_list, str):
        return {}
    # Handles both single/double quotes & spaces
    found = re.findall(r"chain ([A-Z]) and resi (\d+)", contact_list)
    res_chain = {(resi, chain): 1 for chain, resi in found}
    return res_chain

def process_csvs(csv_folders):
    """Aggregate data across all CSVs and folders as needed for plotting."""
    # Dict: base_design -> gene_name -> list of per-model (resi, chain) sets
    all_models = defaultdict(lambda: defaultdict(list))
    for folder in csv_folders:
        for filename in os.listdir(folder):
            if filename.endswith(".csv"):
                df = pd.read_csv(os.path.join(folder, filename))
                df['base_design'] = df['gene_name'].apply(lambda x: x[:-1] if isinstance(x, str) and len(x) > 1 else x)
                for (base_design, gene_name), subdf in df.groupby(['base_design','gene_name']):
                    for contact_list in subdf['contact_list']:
                        res_chain = get_residue_chain_presence(contact_list)
                        all_models[base_design][gene_name].append(res_chain)
    return all_models

def plot_residue_frequency_subplots(all_models, output_folder, max_cols=2):


    os.makedirs(output_folder, exist_ok=True)
    color_palette = [
        'skyblue', 'orange', 'green', 'purple', 'red', 'brown', 'pink', 'olive', 
        'navy', 'teal', 'gold', 'cyan', 'magenta', 'coral', 'gray'
    ]

    for base_design, gene_dict in all_models.items():
        variants = list(gene_dict.keys())
        n_variants = len(variants)
        n_cols = max_cols
        n_rows = math.ceil(n_variants / n_cols)
        # Gather all unique residue numbers and all chain IDs for legend/colors
        all_resi_chain = set()
        for models in gene_dict.values():
            for resi_chain_set in models:
                all_resi_chain.update(resi_chain_set.keys())
        all_residues = sorted({resi for (resi, chain) in all_resi_chain}, key=lambda x: int(x))
        all_chains = sorted({chain for (resi, chain) in all_resi_chain if chain})
        # Assign colors to each chain
        chain_color_dict = {chain: color_palette[i % len(color_palette)] for i, chain in enumerate(all_chains)}

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), squeeze=False)
        fig.suptitle(f'Frequency % by Variant for Base Design {base_design}', fontsize=16, y=0.98)

        for idx, gene_name in enumerate(variants):
            row, col = divmod(idx, n_cols)
            ax = axes[row][col]
            models = gene_dict[gene_name]
            n_models = len(models)
            resi_freq_count = Counter()
            chain_for_resi = dict()

            for m in models:
                for resi_chain in m.keys():
                    resi, chain = resi_chain
                    resi_freq_count[resi] += 1
                    if resi not in chain_for_resi:
                        chain_for_resi[resi] = chain  # assign first chain

            heights = []
            bar_colors = []
            for resi in all_residues:
                freq = resi_freq_count.get(resi, 0) / n_models * 100
                heights.append(freq)
                the_chain = chain_for_resi.get(resi, all_chains[0] if all_chains else 'A')
                bar_colors.append(chain_color_dict.get(the_chain, 'gray'))

            xticks = list(range(len(all_residues)))
            bars = ax.bar(xticks, heights, color=bar_colors, edgecolor='black')
            ax.set_title(gene_name)
            ax.set_xlabel('Residue Position')
            ax.set_ylabel('Frequency %')
            ax.set_xticks(xticks)
            ax.set_xticklabels(all_residues, rotation=45)
            ax.set_ylim(0, 100)

        # Hide unused subplots
        for empty_idx in range(n_variants, n_rows * n_cols):
            row, col = divmod(empty_idx, n_cols)
            axes[row][col].axis('off')

        # Figure-level legend under suptitle
        legend_handles = [plt.Rectangle((0, 0), 1, 1, color=chain_color_dict[c]) for c in all_chains]
        fig.legend(
            legend_handles, all_chains, title='Chain', loc='upper center', 
            bbox_to_anchor=(0.5, 0.92), ncol=min(4, len(all_chains)), frameon=False, fontsize='large'
        )

        plt.tight_layout(rect=[0, 0, 1, 0.88])
        save_path = os.path.join(output_folder, f'{base_design}_residue_frequency_subplots.png')
        plt.savefig(save_path, dpi=300)
        plt.close()


if __name__ == "__main__":
    all_models = process_csvs(csv_folders)
    plot_residue_frequency_subplots(all_models, output_folder, max_cols=2)
    print(f"Subplot figures saved to folder: {output_folder}")
