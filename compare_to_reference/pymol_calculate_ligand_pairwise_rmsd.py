"""
Pairwise RMSD Calculation Script for PyMOL Models

This script calculates and records all pairwise RMSD values within groups of models that share the same gene_name.
Each model is assumed to be named using the format <gene_name>_<ID>, where gene_name may include underscores and the last underscore separates gene_name and the model's unique ID.

For each group with the same gene_name, the script computes the RMSD for every unique model pair and saves the results to a CSV file with the following columns:
    gene_name, model_1, model_2, rmsd

REQUIREMENTS:
- Runs in a PyMOL session with models already loaded and named using the format <gene_name>_<ID>.
- tqdm must be installed: pip install tqdm

USAGE:
    1. Specify a list of reference model names (reference_names) and your desired atom/residue selection substring (selection_substring).
    2. Set the output CSV file path (csv_path)
    3. Load the models and only have the ones you want included in the calculations relative to the references enabled. 
    4. Execute the script:
           run /path/to/ligand_pairwise_rmsd.py
    5. The script will save the CSV file in the specified location with RMSD values for all unique model pairs with the same gene_name.

EXAMPLE OF MAIN EXECUTION SECTION:
    object_names = cmd.get_names("objects",1)
    selection_substring = " and resn PYR"
    csv_path = "./pairwise_rmsd.csv"
    calculate_pairwise_rmsd_and_save(object_names, selection_substring, csv_path)

OUTPUT:
- CSV file with columns: gene_name, model_1, model_2, rmsd
- Each row represents the RMSD between a unique pair of models from the same gene_name group

"""
import csv
import os
from collections import defaultdict
from tqdm import tqdm

# --- User configuration ---
selection_substring = " and resn LIG_B"  # Modify as needed or empty string for whole object
suffix = "_entire_ligand"
csv_path = f"./rmsd_all_to_all_pairwise{suffix}.csv"
# --------------------------

def calculate_pairwise_rmsd_and_save(object_names, substring, csv_filename):
    cmd.remove("hydrogens")
    # Precompute gene_name mapping for each object
    gene_name_dict = {name: name.rsplit('_', 1)[0] for name in object_names}
    # Group models by gene_name
    groups = defaultdict(list)
    for name in object_names:
        gene = gene_name_dict[name]
        groups[gene].append(name)
    results = []
    for gene, models in groups.items():
        n = len(models)
        for i in tqdm(range(n), desc=f"Processing gene {gene}", leave=False):
            for j in range(i+1, n):
                m1 = models[i]
                m2 = models[j]
                cmd.cealign(m2, m1)  # Align structures before RMSD calculation
                try:
                    rmsd = cmd.rms_cur(m1 + substring, m2 + substring)
                except Exception:
                    rmsd = None
                results.append([gene, m1, m2, rmsd])
    # Append mode if file exists, write header only if new
    write_header = not os.path.exists(csv_filename) or os.stat(csv_filename).st_size == 0
    with open(csv_filename, mode='a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(['gene_name', 'model_1', 'model_2', 'rmsd'])
        writer.writerows(results)

# Example main execution - customize before executing in PyMOL
object_names = cmd.get_names("objects",1)
calculate_pairwise_rmsd_and_save(object_names, selection_substring, csv_path)
