"""
RMSD Calculation Script for PyMOL
This script computes the RMSD of protein models relative to a reference model (default: "RelaxClassic") within each gene group, using structures loaded in PyMOL. Outputs a CSV summarizing gene, protocol, model, and RMSD.

REQUIREMENTS:
- PyMOL with models loaded, named as <gene_name>_<ID>.pdb.
- Optional: tqdm for progress bars (pip install tqdm).

USAGE:
- In PyMOL, be in the location of the files you want to analyze.
- (Optional) Set selection_substring to limit RMSD calculation to specific atoms/residues, or leave as "" for all.
- Set csv_path for the output CSV file.
- Execute the script in PyMOL with run /path/to/this_script.py.
- Call the function calculate_rmsd_to_lowest_filtered_pdb(csv_path, selection_substring) to run calculations and write results.

OUTPUT:
- A CSV file listing: gene_name, protocol, model, rmsd for each gene and structural model.
- Each row records the RMSD of a given model compared to the gene's reference ("RelaxClassic") model.
- Models missing the reference are skipped with a warning.

"""
import csv
import os
from collections import defaultdict
import sys

#selection_substring and suffix will be global variables set by piecewise_pymol_analysis.py
# --- User configuration ---
selection_substring = ""  # Modify as needed or empty string for whole object
suffix = ""
csv_path = f"../rmsd_vs_standard{suffix}.csv"
# --------------------------

def find_gene_names():
    gene_names = set()
    files = [f for f in os.listdir('.') if f.endswith('.pdb') and os.path.isfile(f)]
    for file in files:
        name = file.split('_')
        gene_name = name[0]
        if gene_name not in gene_names:
            gene_names.add(gene_name)
    return list(gene_names)

def load_files(gene_name):
    files = [f for f in os.listdir('.') if f.endswith('.pdb') and os.path.isfile(f)]
    filtered = [item for item in files if gene_name in item]
    for file in filtered:
        index = file[-5]
        if 'labmate' in file:
            cmd.load(file, f"{gene_name}_FastRelaxCustom_{index}")
        elif 'std' in file:
            cmd.load(file, f"{gene_name}_FastRelaxStandard_{index}")
        else:
            cmd.load(file, f"{gene_name}_RelaxClassic_{index}")

def find_category(name):
    if 'RelaxClassic' in name:
        return 'RelaxClassic'
    elif 'FastRelaxCustom' in name:
        return 'FastRelaxCustom'
    elif 'FastRelaxStandard' in name:
        return 'FastRelaxStandard'
    else:
        return 'Unknown'

def calculate_rmsd(names,selection_filter): # this is where I left off and need to do more work
    substring = 'RelaxClassic' # Modify as needed
    reference = next((item for item in names if substring in item), None) # get the first one that matches
    others = [item for item in names if item != reference]
    results = []
    for name in others:
        category = find_category(name)
        cmd.align(name,reference, cycles=0, transform=0)
        rmsd_value = cmd.rms_cur(name + selection_filter, reference + selection_filter)
        results.append([category, name, rmsd_value])
    return results

def calculate_rmsd_to_lowest_filtered_pdb(csv_filename,selection_filter):
    compiled_results = []
    genes = find_gene_names()
    for gene in genes:
        load_files(gene)
        cmd.remove("hydrogens")
        all_objects = cmd.get_names("objects", enabled_only=1)
        gene_results = calculate_rmsd(all_objects,selection_filter)
        # Add the gene name to each result
        for result in gene_results:
            result.insert(0, gene)
            compiled_results.append(result)
        cmd.delete('all')
    # Determine file mode: write header only if file does not exist
    write_header = not os.path.exists(csv_filename) or os.stat(csv_filename).st_size == 0
    
    # Write or append to CSV
    with open(csv_filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        if write_header:
            writer.writerow(['gene_name', 'protocol', 'model', 'rmsd'])
        writer.writerows(compiled_results)


# Example main execution - customize before executing in PyMOL
calculate_rmsd_to_lowest_filtered_pdb(csv_path,selection_substring)
